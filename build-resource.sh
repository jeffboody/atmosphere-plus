export RESOURCE=$PWD/app/src/main/assets/resource.bfs

rm $RESOURCE
cd resource

echo BUILD-SHADERS
cd shaders
cd ..

echo ADD-RESOURCES
bfs $RESOURCE blobSet readme.txt

echo CLEANUP
rm shaders/*.spv
cd ..

echo CONTENTS
bfs $RESOURCE blobList
