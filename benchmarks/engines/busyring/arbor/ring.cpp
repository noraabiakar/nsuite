#include <fstream>
#include <iomanip>
#include <iostream>

#include <nlohmann/json.hpp>

#include <arbor/assert_macro.hpp>
#include <arbor/common_types.hpp>
#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/simple_sampler.hpp>
#include <arbor/simulation.hpp>
#include <arbor/recipe.hpp>
#include <arbor/version.hpp>

#include <arborenv/concurrency.hpp>
#include <arborenv/gpu_env.hpp>

#include "parameters.hpp"

#ifdef ARB_MPI_ENABLED
#include <mpi.h>
#include <arborenv/with_mpi.hpp>
#endif

using arb::cell_gid_type;
using arb::cell_lid_type;
using arb::cell_size_type;
using arb::cell_member_type;
using arb::cell_kind;
using arb::time_type;
using arb::cell_probe_address;

// Writes voltage trace as a json file.
void write_trace_json(std::string fname, const arb::trace_data<double>& trace);

// Generate a cell.
arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params);

class ring_recipe: public arb::recipe {
public:
    ring_recipe(ring_params params):
        num_cells_(params.num_cells),
        min_delay_(params.min_delay),
        params_(params)
    {
        gprop.default_parameters = arb::neuron_parameter_defaults;
        gprop.default_parameters.temperature_K = 35 + 273.15;
        gprop.default_parameters.init_membrane_potential = -65;

        gprop.default_parameters.reversal_potential_method["k"] = "nernst/k";
        gprop.default_parameters.reversal_potential_method["na"] = "nernst/na";
    }

    cell_size_type num_cells() const override {
        return num_cells_;
    }

    arb::util::unique_any get_cell_description(cell_gid_type gid) const override {
        return branch_cell(gid, params_.cell);
    }

    cell_kind get_cell_kind(cell_gid_type gid) const override {
        return cell_kind::cable;
    }

    // Each cell has one spike detector (at the soma).
    cell_size_type num_sources(cell_gid_type gid) const override {
        return 1;
    }

    // The cell has one target synapse, which will be connected to cell gid-1.
    cell_size_type num_targets(cell_gid_type gid) const override {
        return params_.cell.synapses;
    }

    arb::util::any get_global_properties(cell_kind kind) const override {
        return gprop;
    }

    // Each cell has one incoming connection, from cell with gid-1,
    // and fan_in-1 random connections with very low weight.
    std::vector<arb::cell_connection> connections_on(cell_gid_type gid) const override {
        std::vector<arb::cell_connection> cons;
        const auto ncons = params_.cell.synapses;
        cons.reserve(ncons);

        const auto s = params_.ring_size;
        const auto group = gid/s;
        const auto group_start = s*group;
        const auto group_end = std::min(group_start+s, num_cells_);
        cell_gid_type src = gid==group_start? group_end-1: gid-1;
        cons.push_back(arb::cell_connection({src, 0}, {gid, 0}, event_weight_, min_delay_));

        // Used to pick source cell for a connection.
        std::uniform_int_distribution<cell_gid_type> dist(0, num_cells_-2);
        // Used to pick delay for a connection.
        std::uniform_real_distribution<float> delay_dist(0, 2*min_delay_);
        auto src_gen = std::mt19937(gid);
        for (unsigned i=1; i<ncons; ++i) {
            // Make a connection with weight 0.
            // The source is randomly picked, with no self connections.
            src = dist(src_gen);
            if (src==gid) ++src;
            const float delay = min_delay_+delay_dist(src_gen);
            //const float delay = min_delay_;
            cons.push_back(
                arb::cell_connection({src, 0}, {gid, i}, 0.f, delay));
            std::cout << src << " -> " << gid << std::endl;
        }

        return cons;
    }

    // Return one event generator on the first cell of each ring.
    // This generates a single event that will kick start the spiking on the sub-ring.
    std::vector<arb::event_generator> event_generators(cell_gid_type gid) const override {
        std::vector<arb::event_generator> gens;
        if (gid%params_.ring_size == 0) {
            gens.push_back(
                arb::explicit_generator(
                    arb::pse_vector{{{gid, 0}, 1.0, event_weight_}}));
        }
        return gens;
    }

    // There is one probe (for measuring voltage at the soma) on the cell.
    cell_size_type num_probes(cell_gid_type gid)  const override {
        return 1;
    }

    arb::probe_info get_probe(cell_member_type id) const override {
        // Get the appropriate kind for measuring voltage.
        cell_probe_address::probe_kind kind = cell_probe_address::membrane_voltage;
        // Measure at the soma.
        arb::mlocation loc{0, 0.0};

        return arb::probe_info{id, kind, cell_probe_address{loc, kind}};
    }

    void add_ion(const std::string& ion_name, int charge, double init_iconc, double init_econc, double init_revpot) {
        gprop.add_ion(ion_name, charge, init_iconc, init_econc, init_revpot);
    }

private:
    cell_size_type num_cells_;
    double min_delay_;
    ring_params params_;
    arb::cable_cell_global_properties gprop;

    float event_weight_ = 0;
};

struct cell_stats {
    using size_type = unsigned;
    size_type ncells = 0;
    size_type nbranch = 0;
    size_type ncomp = 0;

    cell_stats(arb::recipe& r) {
#ifdef ARB_MPI_ENABLED
        int nranks, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nranks);
        ncells = r.num_cells();
        size_type cells_per_rank = ncells/nranks;
        size_type b = rank*cells_per_rank;
        size_type e = (rank==nranks-1)? ncells: (rank+1)*cells_per_rank;
        size_type nbranch_tmp = 0;
        size_type ncomp_tmp = 0;
        for (size_type i=b; i<e; ++i) {
            auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
            nbranch_tmp += c.num_segments();
            ncomp_tmp += c.num_compartments();
        }
        MPI_Allreduce(&nbranch_tmp, &nbranch, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&ncomp_tmp, &ncomp, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
        ncells = r.num_cells();
        for (size_type i=0; i<ncells; ++i) {
            auto c = arb::util::any_cast<arb::cable_cell>(r.get_cell_description(i));
            nbranch += c.num_segments();
            ncomp += c.num_compartments();
        }
#endif
    }

    friend std::ostream& operator<<(std::ostream& o, const cell_stats& s) {
        return o << "cell stats: "
                 << s.ncells << " cells; "
                 << s.nbranch << " branches; "
                 << s.ncomp << " compartments.";
    }
};

int main(int argc, char** argv) {
    try {
        bool root = true;

        auto params = read_options(argc, argv);

        arb::proc_allocation resources;
        if (auto nt = arbenv::get_env_num_threads()) {
            resources.num_threads = nt;
        }
        else {
            resources.num_threads = arbenv::thread_concurrency();
        }

#ifdef ARB_MPI_ENABLED
        arbenv::with_mpi guard(argc, argv, false);
        resources.gpu_id = arbenv::find_private_gpu(MPI_COMM_WORLD);
        auto context = arb::make_context(resources, MPI_COMM_WORLD);
        root = arb::rank(context) == 0;
#else
        resources.gpu_id = arbenv::default_gpu();
        auto context = arb::make_context(resources);
#endif

#ifdef ARB_PROFILE_ENABLED
        arb::profile::profiler_initialize(context);
#endif

        // Print a banner with information about hardware configuration
        if (root) {
            std::cout << "gpu:      " << (has_gpu(context)? "yes": "no") << "\n";
            std::cout << "threads:  " << num_threads(context) << "\n";
            std::cout << "mpi:      " << (has_mpi(context)? "yes": "no") << "\n";
            std::cout << "ranks:    " << num_ranks(context) << "\n" << std::endl;
        }

        arb::profile::meter_manager meters;
        meters.start(context);

        // Create an instance of our recipe.
        ring_recipe recipe(params);
        recipe.add_ion("h", 1, 1.0, 1.0, -34.4);
        recipe.add_ion("no", 1, 1.0, 1.0, 0);

        cell_stats stats(recipe);
        if (root) std::cout << stats << "\n";

        //arb::partition_hint_map hints;
        //hints[cell_kind::cable1d_neuron].cpu_group_size = 4;
        //auto decomp = arb::partition_load_balance(recipe, context, hints);
        auto decomp = arb::partition_load_balance(recipe, context);

        // Construct the model.
        arb::simulation sim(recipe, decomp, context);

        // Set up the probe that will measure voltage in the cell.

        // This is where the voltage samples will be stored as (time, value) pairs
        arb::trace_data<double> voltage;
        if (params.record_voltage) {
            // The id of the only probe on the cell:
            // the cell_member type points to (cell 0, probe 0)
            auto probe_id = cell_member_type{0, 0};
            // The schedule for sampling is 10 samples every 1 ms.
            auto sched = arb::regular_schedule(params.dt);
            // Now attach the sampler at probe_id, with sampling schedule sched, writing to voltage
            sim.add_sampler(arb::one_probe(probe_id), sched, arb::make_simple_sampler(voltage));
        }

        // Set up recording of spikes to a vector on the root process.
        std::vector<arb::spike> recorded_spikes;
        if (root) {
            sim.set_global_spike_callback(
                [&recorded_spikes](const std::vector<arb::spike>& spikes) {
                    recorded_spikes.insert(recorded_spikes.end(), spikes.begin(), spikes.end());
                });
        }

        meters.checkpoint("model-init", context);

        // Run the simulation.
        if (root) std::cout << "running simulation" << std::endl;
        sim.set_binning_policy(arb::binning_kind::regular, params.dt);
        sim.run(params.duration, params.dt);

        meters.checkpoint("model-run", context);

        auto ns = sim.num_spikes();

        // Write spikes to file
        if (root) {
            std::cout << "\n" << ns << " spikes generated at rate of "
                      << params.duration/ns << " ms between spikes\n";
            std::ofstream fid(params.odir + "/" + params.name + "_spikes.gdf");
            if (!fid.good()) {
                std::cerr << "Warning: unable to open file spikes.gdf for spike output\n";
            }
            else {
                char linebuf[45];
                for (auto spike: recorded_spikes) {
                    auto n = std::snprintf(
                        linebuf, sizeof(linebuf), "%u %.4f\n",
                        unsigned{spike.source.gid}, float(spike.time));
                    fid.write(linebuf, n);
                }
            }
        }

        // Write the samples to a json file samples were stored on this rank.
        if (voltage.size()>0u) {
            std::string fname = params.odir + "/" + params.name + "_voltages.json";
            write_trace_json(fname, voltage);
        }

        auto report = arb::profile::make_meter_report(meters, context);
        if (root) std::cout << report;
    }
    catch (std::exception& e) {
        std::cerr << "exception caught in ring miniapp: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

void write_trace_json(std::string fname, const arb::trace_data<double>& trace) {
    nlohmann::json json;
    json["name"] = "ring demo";
    json["units"] = "mV";
    json["cell"] = "0.0";
    json["probe"] = "0";

    auto& jt = json["data"]["time"];
    auto& jy = json["data"]["voltage"];

    for (const auto& sample: trace) {
        jt.push_back(sample.t);
        jy.push_back(sample.v);
    }

    std::ofstream file(fname);
    file << std::setw(1) << json << "\n";
}

// Helper used to interpolate in branch_cell.
template <typename T>
double interp(const std::array<T,2>& r, unsigned i, unsigned n) {
    double p = i * 1./(n-1);
    double r0 = r[0];
    double r1 = r[1];
    return r[0] + p*(r1-r0);
}

arb::sample_tree read_swc(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("unable to open SWC file: "+path);

    return arb::swc_as_sample_tree(arb::parse_swc_file(f));
}

arb::cable_cell branch_cell(arb::cell_gid_type gid, const cell_parameters& params) {
    auto tree = read_swc("/users/abiakarn/git/nsuite_fork/benchmarks/engines/busyring/arbor/purkinje_mouse_original.swc");

    arb::label_dict dict;
    using arb::reg::tagged;
    dict.set("soma", tagged(1));
    dict.set("axon", tagged(2));

    std::vector<arb::msample> axon_samples = {
            {{0, 0, 0,   0.5}, 2},
            {{0, 0, 17,  0.5}, 2},
            {{0, 0, 21,  0.5}, 2},
            {{0, 0, 121, 0.5}, 2},
            {{0, 0, 125, 0.5}, 2},
            {{0, 0, 225, 0.5}, 2},
            {{0, 0, 229, 0.5}, 2},
            {{0, 0, 329, 0.5}, 2},
            {{0, 0, 333, 0.5}, 2},
            {{0, 0, 433, 0.5}, 2}
    };
    tree.append(0, axon_samples);
    arb::morphology morpho(tree);

    arb::cable_cell c = make_cable_cell(morpho, dict, false);

    // Paint mechanisms on the soma
    {
        arb::mechanism_desc Leak("Leak");
        Leak["e"] = -61;
        Leak["gmax"] = 0.001;

        arb::mechanism_desc Nav1_6("Nav1_6");
        Nav1_6["gbar"] = 0.18596749324385001;

        arb::mechanism_desc Kv1_1("Kv1_1");
        Kv1_1["gbar"] = 0.0029172006269699998;

        arb::mechanism_desc Kv3_4("Kv3_4");
        Kv3_4["gkbar"] = 0.069972751903779995;

        arb::mechanism_desc Kir2_3("Kir2_3");
        Kir2_3["gkbar"] = 2.322613156e-05;

        arb::mechanism_desc Kca1_1("Kca1_1");
        Kca1_1["gbar"] = 0.01197387128516;

        arb::mechanism_desc Kca2_2("Kca2_2");
        Kca2_2["gkbar"] = 0.0013377920303699999;

        arb::mechanism_desc Kca3_1("Kca3_1");
        Kca3_1["gkbar"] = 0.01388910824701;

        arb::mechanism_desc Cav2_1("Cav2_1");
        Cav2_1["pcabar"] = 0.00020306777733000001;

        arb::mechanism_desc Cav3_1("Cav3_1");
        Cav3_1["pcabar"] = 5.1352684600000001e-06;

        arb::mechanism_desc Cav3_2("Cav3_2");
        Cav3_2["gcabar"] = 0.00070742370991999995;

        arb::mechanism_desc Cav3_3("Cav3_3");
        Cav3_3["pcabar"] = 0.00034648446559000002;

        arb::mechanism_desc HCN1("HCN1");
        HCN1["gbar"] = 0.0016391306742300001;

        arb::mechanism_desc cdp5("cdp5");
        cdp5["TotalPump"] = 2e-08;

        arb::mechanism_desc kv15("Kv1_5");
        kv15["gKur"] = 0.00011449636712999999;

        arb::mechanism_desc kv33("Kv3_3");
        kv33["gbar"] = 0.01054618632087;

        arb::mechanism_desc kv43("Kv4_3");
        kv43["gkbar"] = 0.00147529033238;

        c.paint("soma", Leak);
        c.paint("soma", Nav1_6);
        c.paint("soma", Kv1_1);
        c.paint("soma", Kv3_4);
        c.paint("soma", Kir2_3);
        c.paint("soma", Kca1_1);
        c.paint("soma", Kca2_2);
        c.paint("soma", Kca3_1);
        c.paint("soma", Cav2_1);
        c.paint("soma", Cav3_1);
        c.paint("soma", Cav3_2);
        c.paint("soma", Cav3_3);
        c.paint("soma", HCN1);
        c.paint("soma", cdp5);

        c.paint("soma", kv15);
        c.paint("soma", kv33);
        c.paint("soma", kv43);
    }

    // Paint mechanisms on the axon
    {
        arb::mechanism_desc Leak("Leak");
        Leak["e"] = -61;
        Leak["gmax"] = 0.00029999999999999997;

        arb::mechanism_desc Nav1_6("Nav1_6");
        Nav1_6["gbar"] = 0.027565014930900002;

        arb::mechanism_desc Kv3_4("Kv3_4");
        Kv3_4["gkbar"] = 0.029949599085;

        arb::mechanism_desc Cav2_1("Cav2_1");
        Cav2_1["pcabar"] = 0.00026695621961;

        arb::mechanism_desc Cav3_1("Cav3_1");
        Cav3_1["pcabar"] = 1.487097008e-05;

        arb::mechanism_desc cdp5("cdp5");
        cdp5["TotalPump"] = 4.9999999999999998e-07;

        c.paint("axon", Leak);
        c.paint("axon", Nav1_6);
        c.paint("axon", Kv3_4);
        c.paint("axon", Cav2_1);
        c.paint("axon", Cav3_1);
        c.paint("axon", cdp5);
    }

    // Set ion parameters, ra and cm on soma
    arb::soma_segment* soma = c.soma();
    soma->parameters.axial_resistivity = 122;
    soma->parameters.membrane_capacitance = 0.0077000000000000002;

//        soma->parameters.ion_data["k"].init_reversal_potential = -88; // not being read in Neuron
    soma->parameters.ion_data["k"].init_int_concentration = 54.4;
    soma->parameters.ion_data["k"].init_ext_concentration = 2.5;

//        soma->parameters.ion_data["na"].init_reversal_potential = 80; // not being read in Neuron
    soma->parameters.ion_data["na"].init_int_concentration = 10;
    soma->parameters.ion_data["na"].init_ext_concentration = 140;

//        soma->parameters.ion_data["ca"].init_reversal_potential = 137.5; // not used anywhere
    soma->parameters.ion_data["ca"].init_int_concentration = 0.00005000;
    soma->parameters.ion_data["ca"].init_ext_concentration = 2.0;

    // Set ion parameters, ra and cm on axon
    arb::cable_segment* axon = c.cable(c.num_segments()-1);
    axon->parameters.axial_resistivity = 122;
    axon->parameters.membrane_capacitance = 0.0077000000000000002;

    axon->parameters.ion_data["k"].init_int_concentration = 54.4;
    axon->parameters.ion_data["k"].init_ext_concentration = 2.5;

    axon->parameters.ion_data["na"].init_int_concentration = 10;
    axon->parameters.ion_data["na"].init_ext_concentration = 140;

    axon->parameters.ion_data["ca"].init_int_concentration = 0.00005000;
    axon->parameters.ion_data["ca"].init_ext_concentration = 2.0;

    // Discretize dendrites: 1 compartment per branch; axon: 9 compartments
    for (std::size_t i=1; i<c.num_segments(); ++i) {
        arb::cable_segment* branch = c.cable(i);
        branch->set_compartments(10);
    }
    c.cable(c.num_segments()-1)->set_compartments(90);

    // Add spike threshold detector at the soma.
    c.place(arb::mlocation{0,0}, arb::threshold_detector{10});

    // Add a synapse to the mid point of the first dendrite.
    c.place(arb::mlocation{1, 0.5}, "expsyn");

    // Add additional synapses that will not be connected to anything.
    for (unsigned i=1u; i<params.synapses; ++i) {
        c.place(arb::mlocation{1, 0.5}, "expsyn");
    }

    return c;
}
