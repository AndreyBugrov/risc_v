uname -r > device_info_$1.txt
lsb_release -a >> device_info_$1.txt
cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies >> device_info_$1.txt
lscpu | grep 'CPU(s):' >> device_info_$1.txt
cat /proc/cpuinfo | grep 'rv' | tail -1 >> device_info_$1.txt
cat /proc/cpuinfo | grep 'isa-ext' | tail -1 >> device_info_$1.txt
getconf -a | grep CACHE >> device_info_$1.txt
free -m >> device_info_$1.txt
g++ --version | head -1 >> device_info_$1.txt
