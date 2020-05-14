cp /usr/bin/gcc /usr/bin/gcc-orig
cp /usr/bin/g++ /usr/bin/g++-orig
cp /usr/bin/ld /usr/bin/ld-orig
cp /usr/bin/cpp /usr/bin/cpp-orig
cp /usr/bin/objcopy /usr/bin/objcopy-orig

rm /usr/bin/gcc /usr/bin/ld /usr/bin/g++ /usr/bin/cpp /usr/bin/objcopy

ln -s /opt/bolos-env/gcc-arm-none-eabi/bin/arm-none-eabi-gcc /usr/bin/gcc
ln -s /opt/bolos-env/gcc-arm-none-eabi/bin/arm-none-eabi-g++ /usr/bin/g++
ln -s /opt/bolos-env/gcc-arm-none-eabi/bin/arm-none-eabi-ld /usr/bin/ld
ln -s /opt/bolos-env/gcc-arm-none-eabi/bin/arm-none-eabi-cpp /usr/bin/cpp
ln -s /opt/bolos-env/gcc-arm-none-eabi/bin/arm-none-eabi-objcopy /usr/bin/objcopy
