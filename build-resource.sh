export RESOURCE=$PWD/app/src/main/assets/resource.bfs

rm $RESOURCE
cd resource

echo BUILD-SHADERS
cd shaders
glslangValidator -V planet.frag -o planet_frag.spv
glslangValidator -V planet.vert -o planet_vert.spv
glslangValidator -V sky_flat.frag -o sky_flat_frag.spv
glslangValidator -V sky_flat.vert -o sky_flat_vert.spv
glslangValidator -V sky_atmo.frag -o sky_atmo_frag.spv
glslangValidator -V sky_atmo.vert -o sky_atmo_vert.spv
cd ..

echo ADD-RESOURCES
bfs $RESOURCE blobSet readme.txt
bfs $RESOURCE blobSet models/Sphere.glb
bfs $RESOURCE blobSet shaders/planet_frag.spv
bfs $RESOURCE blobSet shaders/planet_vert.spv
bfs $RESOURCE blobSet shaders/sky_flat_frag.spv
bfs $RESOURCE blobSet shaders/sky_flat_vert.spv
bfs $RESOURCE blobSet shaders/sky_atmo_frag.spv
bfs $RESOURCE blobSet shaders/sky_atmo_vert.spv

echo CLEANUP
rm shaders/*.spv
cd ..

echo ADD-VKK
cd app/src/main/cpp/libvkk/ui/resource
./build-resource.sh $RESOURCE
cd ../../../../../../..

echo CONTENTS
bfs $RESOURCE blobList
