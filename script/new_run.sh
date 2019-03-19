# mkdir analysis & cd analysis
# ../script/new_run.sh <test_name> ..

mkdir $1
cd $1
basedir=$2
ln -s ${basedir}/../script .
ln -s ${basedir}/../data . 
ln -s ${basedir}/../resources .
ln -s ${basedir}/../main.nf .
ln -s ${basedir}/../nextflow.config .
# ln -s ${basedir}/../result .
