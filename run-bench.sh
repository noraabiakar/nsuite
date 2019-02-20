usage() {
    echo
    echo "nsuite benchmark runner:"
    echo
    echo "run with one, and only one, of the following options:"
    echo "   arbor  : run arbor benchmarks"
    echo "   neuron : run neuron benchmarks"
    echo "   coreneuron : run coreneuron benchmarks"
    echo

    exit 1;
}

run_arb=false
run_nrn=false
run_corenrn=false

models="ring"
configs="small"

# Load some utility functions.
source ./scripts/environment.sh
source ./scripts/util.sh
default_environment

# parse arguments
while [ "$1" != "" ]
do
    case $1 in
        arbor )
            run_arb=true
            ;;

        neuron )
            run_nrn=true
            ;;

        coreneuron )
            run_corenrn=true
            ;;
        --model )
            shift
            models="$1"
            ;;
        --config )
            shift
            configs="$1"
            ;;

        * )
            echo "unknown option '$1'"
            usage
    esac
    shift
done

if [ "$run_arb" == "true" ]; then
    [ ! -f "$ns_base_path/config/env_arbor.sh" ] &&  err "Arbor must be installed to run Arbor benchmarks." && run_arb=false
fi
if [ "$run_nrn" == "true" ]; then
    [ ! -f "$ns_base_path/config/env_neuron.sh" ] &&  err "NEURON must be installed to run NEURON benchmarks." && run_nrn=false
fi
if [ "$run_corenrn" == "true" ]; then
    [ ! -f "$ns_base_path/config/env_coreneuron.sh" ] &&  err "CoreNeuron must be installed to run CoreNeuron benchmarks." && run_corenrn=false
fi

# TODO: this has to go into the configuration environment setup scripts
export ARB_NUM_THREADS=$[ $ns_threads_per_core * $ns_cores_per_socket ]

echo "---- Platform ----"
echo "configuration:     $ns_environment"
echo "platform:          $ns_system ($(uname -or))"
echo "cores per socket:  $ns_cores_per_socket"
echo "threads per core:  $ns_threads_per_core"
echo "threads:           $ARB_NUM_THREADS"
echo "sockets:           $ns_sockets"
echo "mpi:               $ns_with_mpi"
echo

echo "---- Application ----"
echo "Arbor:      $run_arb"
echo "NEURON:     $run_nrn"
echo "CoreNeuron: $run_corenrn"
echo

echo "---- Benchmarks ----"
echo

mkdir -p "$ns_input_path"
for model in $models
do
    model_config_path="$ns_base_path/benchmarks/models/$model"

    cd "$model_config_path"

    for config in $configs
    do

        echo $model-$config
        echo

        model_input_path="$ns_input_path/benchmarks/$model/$config"
        model_output_path="$ns_output_path/benchmarks/$model/$config"

        ./config.sh $config "$ns_base_path" "$model_input_path" "$model_output_path" "$ns_base_path/config"

        # todo: hoist check for env file outside loop, which would unset any simulation engine that has not been installed
        if [ "$run_arb" == "true" ]; then
            echo "  arbor"
            "$model_input_path/run_arb.sh"
        fi
        if [ "$run_nrn" == "true" ]; then
            echo "  neuron"
            "$model_input_path/run_nrn.sh"
        fi
        if [ "$run_corenrn" == "true" ]; then
            echo "  coreneuron"
            "$model_input_path/run_corenrn.sh"
        fi
    done
done
