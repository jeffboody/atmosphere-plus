export RESOURCE=$PWD/app/src/main/assets/resource.bfs

rm $RESOURCE
cd resource

echo BUILD-SHADERS
cd shaders
glslangValidator -V planet.frag -o planet_frag.spv
glslangValidator -V planet.vert -o planet_vert.spv
cd ..

echo ADD-RESOURCES
bfs $RESOURCE blobSet readme.txt
bfs $RESOURCE blobSet models/Sphere.glb
bfs $RESOURCE blobSet shaders/planet_frag.spv
bfs $RESOURCE blobSet shaders/planet_vert.spv

echo CLEANUP
rm shaders/*.spv
cd ..

echo CONTENTS
bfs $RESOURCE blobList
