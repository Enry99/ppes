ppes_path=$(pwd)

if !(grep -q ${ppes_path} ~/.bashrc)
then
echo "export PATH=${ppes_path}/ppes:\${PATH}" >> ~/.bashrc
echo "export PYTHONPATH=${ppes_path}:\${PYTHONPATH}" >> ~/.bashrc
else
echo "ppes already found on PATH. Please remove it from .bashrc if you want to change the script location."
fi

chmod +x ppes/ppes


source ~/.bashrc
