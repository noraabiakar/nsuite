SCRIPT_DIR="$(dirname "$(realpath "$0")")"

if  [ "$1" != "" ]; then
  PLATFORM=$1 
else
  PLATFORM="mc"
fi

echo "Building Arbor benchmarks"
jobid=`sbatch --wait $SCRIPT_DIR/batch_build_${PLATFORM}.sh`
ret=$?

echo -n "Result: "
if [ $ret -ne 0 ]; then
  echo "BUILD FAILED."
  echo "============"
  echo "   OUTPUT"
  echo "============"
  cat build_${PLATFORM}.err
  exit $ret
else 
  echo "BUILD SUCCESSFUL."
  echo "============"
  echo "   OUTPUT"
  echo "============"
  cat build_${PLATFORM}.out
fi

echo ""
echo "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
echo ""
echo "Running Arbor benchmarks"
jobid=`sbatch --wait $SCRIPT_DIR/batch_run_${PLATFORM}.sh`
ret=$?

echo -n "Result: "
if [ $ret -ne 0 ]; then
  echo "RUNNING BENCHMARKS FAILED."
  echo "============"
  echo "   OUTPUT"
  echo "============"
  cat run_${PLATFORM}.err
  exit $ret
else 
  echo "RUNNING BENCHMARKS SUCCESSFUL."
  echo "============"
  echo "   OUTPUT"
  echo "============"
  cat run_${PLATFORM}.out
fi
