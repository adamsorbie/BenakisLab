echo "Enter folder name: " 
read folder 

mkdir -p $folder 
cd $folder 
curl -LOk https://github.com/Lagkouvardos/Rhea/archive/master.zip 
unzip master.zip
rm -f master.zip 