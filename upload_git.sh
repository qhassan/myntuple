cp  ~/public/ttree/CMSSW_7_0_6_patch1/src/Test/MiniAnalyzer/plugins/* plugins/
cd plugins
git add .
cd ..
git commit -m "$1"
git push origin master 
