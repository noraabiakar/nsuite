usage() {
    echo
    echo "nsuite installer options:"
    echo
    echo "   arbor       : build Arbor"
    echo "   neuron      : build NEURON"
    echo "   coreneuron  : build CoreNEURON"
    echo "   all         : build all simulators"
    echo "   -e filename : source filename before building"
    echo
    echo "examples:"
    echo
    echo "install only Arbor:"
    echo "$ install arbor"
    echo
    echo "install Arbor, NEURON and CoreNEURON:"
    echo "$ install all"
    echo
    echo "install NEURON using environment configured in config.sh:"
    echo "$ install neuron -e config.sh"
    echo
}

# Load some utility functions.
source ./scripts/util.sh
source ./scripts/environment.sh

# Set up default environment variables
default_environment

# parse arguments
while [ "$1" != "" ]
do
    case $1 in
        arbor )
            ns_build_arbor=true
            ;;
        neuron )
            ns_build_neuron=true
            ;;
        coreneuron )
            ns_build_coreneuron=true
            ;;
        all )
            ns_build_arbor=true
            ns_build_neuron=true
            ns_build_coreneuron=true
            ;;
        -e )
            shift
            ns_environment=$1
            ;;
        * )
            echo "unknown option '$1'"
            usage
            exit 1
    esac
    shift
done

# Run a user supplied configuration script if it was provided with the -e flag.
# This will make changes to the configuration variables ns_* set in environment()
if [ "$ns_environment" != "" ]; then
    echo "using additional configuration: $ns_environment"
    if [ ! -f "$ns_environment" ]; then
        err "file '$ns_environment' not found"
        exit 1
    fi
    source "$ns_environment"
    echo
fi

echo "---- TARGETS ----"
echo "build arbor:       $ns_build_arbor"
echo "build neuron:      $ns_build_neuron"
echo "build coreneuron:  $ns_build_coreneuron"
echo
echo "---- PATHS ----"
echo "working path:  $ns_base_path"
echo "install path:  $ns_install_path"
echo "build path:    $ns_build_path"
echo "input path:    $ns_input_path"
echo "output path:   $ns_output_path"
echo
echo "---- SYSTEM ----"
echo "system:        $ns_system"
echo "using mpi:     $ns_with_mpi"
echo "C compiler:    $ns_cc"
echo "C++ compiler:  $ns_cxx"
echo "python:        $ns_python"
echo
echo "---- ARBOR ----"
echo "repo:          $ns_arb_repo"
echo "branch:        $ns_arb_branch"
echo "arch:          $ns_arb_arch"
echo "gpu:           $ns_arb_with_gpu"
echo "vectorize:     $ns_arb_vectorize"
echo
echo "---- NEURON ----"
echo "tarball:       $ns_nrn_tarball"
echo "url:           $ns_nrn_url"
echo "repo:          $ns_nrn_git_repo"
echo "branch:        $ns_nrn_branch"
echo
echo "---- CoreNEURON ----"
echo "repo:          $ns_cnrn_git_repo"
echo "sha:           $ns_cnrn_sha"

mkdir -p "$ns_build_path"

export CC="$ns_cc"
export CXX="$ns_cxx"

[ "$ns_build_arbor"  = true ] && echo && source "$ns_base_path/scripts/build_arbor.sh"
cd "$ns_base_path"
[ "$ns_build_neuron" = true ] && echo && source "$ns_base_path/scripts/build_neuron.sh"
cd "$ns_base_path"
[ "$ns_build_coreneuron" = true ] && echo && source "$ns_base_path/scripts/build_coreneuron.sh"
cd "$ns_base_path"

echo
echo "Installation finished"
echo

# Find and record the python and binary paths.
find_paths python_path site-packages
find_paths bin_path bin

echo "python paths: $python_path"
echo "bin paths:    $bin_path"

#config_path="${ns_base_path}/config"
#config_file="${config_path}/env.sh"
#mkdir -p "$config_path"

#echo "export PATH=\"${ns_install_path}/bin:\${PATH}\""  > "$config_file"
#echo "export PYTHONPATH=\"${ns_base_path}/common/python:\${PYTHONPATH}\"" >> "$config_file"
#echo "export PYTHONPATH=\"$python_path\$PYTHONPATH\""   >> "$config_file"
#echo "export PATH=\"$bin_path\$PATH\""                  >> "$config_file"
#echo "source \"$ns_base_path/scripts/environment.sh\""  >> "$config_file"
#echo "default_environment"                              >> "$config_file"
#if [ "$ns_environment" != "" ]; then
#    full_env=$(full_path "$ns_environment")
#    echo "source \"$full_env\""                         >> "$config_file"
#fi
