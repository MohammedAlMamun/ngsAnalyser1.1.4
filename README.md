# shiny-electron-ngsAnalysis
Electron-Shiny standalone app to bring DNA sequencing data analysis to the point of your mouse-clicks!
## prerequisites
macOS 10.13 (High Sierra) and higher<br/>
R 4.1.1 (https://cran.r-project.org/bin/macosx/)<br/>
Python 3.7.7 (https://www.python.org/downloads/release/python-377/)<br/>
XQuartz latest (https://www.xquartz.org/)<br/>
<br/>
*you may use the links above to download appropriate versions; double click the .pkg file and follow the on screen instructions to install. XQuartz .pkg is inside the .dmg file.
<br/>
## install pyhton virtual environment
### open terminal window
> click top menu's Go button and choose utilities<br/>
> inside utilities double-click Terminal<br/>
> it will open up a Terminal window
### check python version
> paste this line below in the terminal prompt<br/>
python --version<br/>
<br/>
> this will give you the following result<br/>
Python 2.7.16<br/>
<br/>
> this is the mac's system python<br/>
> but we will use python3 and that's why we installed Python 3.7.7, now we will check if the installation was fine<br/>
<br/>
> copy and paste the following command in your terminal prompt<br/>
python3 --version<br/>
<br/>
> if everything okay, you will have the following resul<br/>
Python 3.7.7<br/>

### install virtual environment

> Now paste the following command in your terminal prompt<br/>
pip3 install virtualenv<br/>
<br/>
> this will install virtual environment and you can confirm the installation by pasting the following command<br/>
virtualenv<br/>
<br/>
> this will produce the following<br/>
<br/>
usage: virtualenv [--version] [--with-traceback] [-v | -q] ... <br/> 
> and the system will exit with error message<br/>
virtualenv: error: the following arguments are required: dest<br/>
<br/>
<br/>
> that's it. we are now ready to launch ngsAnalyser.






