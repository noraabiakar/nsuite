SCRIPT_DIR="$(dirname "$(realpath "$0")")"

jobid=`timeout -s 9 1m sbatch --parsable $SCRIPT_DIR/batch_run.sh`
ret=$?
echo "Submitted job $jobid"

if [ $ret -ne 0 ]; then
  echo "Sbatch failed."
  exit $ret
fi

status=""

# sacct output may be empty until the job shows up.
while [ "$status" == "PENDING" -o "$status" == "RUNNING" -o "$status" == "" ]; do
  status=`sacct -j ${jobid} -o State -n -P | head -n 1`
  echo "Status $status"
  sleep 5
done

echo ""

echo "----- Result: -----"
if [ "$status" == "COMPLETED" ]; then
  echo "Test passed."
else
  echo "Test FAILED. (Status: `sacct -j ${jobid} -oState,ExitCode -n | head -n 1`)"
fi

echo ""

if [ "$status" != "COMPLETED" ]; then
  exit 2
fi

