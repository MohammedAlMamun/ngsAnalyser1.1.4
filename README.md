# shiny-electron-ngsAnalysis
Electron-Shiny standalone app to bring DNA sequencing data analysis to the point of your mouse-clicks!
## prerequisites
macOS 10.13 (High Sierra) and higher<br/>
R 4.1.1 (https://cran.r-project.org/bin/macosx/)<br/>
Python 3.7.7 (https://www.python.org/downloads/release/python-377/)<br/>
XQuartz latest (https://www.xquartz.org/)<br/>
<br/>
_you may use the links above to download appropriate versions; double click the .pkg file and follow the on screen instructions to install. XQuartz .pkg is inside the .dmg file._
<br/>
## install pyhton virtual environment
### open terminal window
> click top menu's Go button and choose utilities<br/>
> inside utilities double-click Terminal<br/>
> it will open up a Terminal window
### check python version
> paste this line below in the terminal prompt<br/>
`python --version`<br/>

> this will give you the following result<br/>
`Python 2.7.16`<br/>

> this is the mac's system python<br/>
> but we will use python3 and that's why we installed Python 3.7.7, now we will check if the installation was fine<br/>

> copy and paste the following command in your terminal prompt<br/>
`python3 --version`<br/>

> if everything okay, you will have the following resul<br/>
`Python 3.7.7`<br/>

### install virtual environment

> Now paste the following command in your terminal prompt<br/>
`pip3 install virtualenv`<br/>

> this will install virtual environment and you can confirm the installation by pasting the following command<br/>
`virtualenv`<br/>

> this will produce the following<br/>
`usage: virtualenv [--version] [--with-traceback] [-v | -q] ...` <br/> 

> and the system will exit with error message<br/>
`virtualenv: error: the following arguments are required: dest`<br/>

> that's it. we are now ready to launch ngsAnalyser.

## install ngsAnalyser
Download or get the shared <br/>
ngsAnalyser-1.1.4.dmg <br/>
double click the file and copy the app file to your Applications folder.<br/>

_now you can simply double-click the app icon and it should open a new window on your desktop._








