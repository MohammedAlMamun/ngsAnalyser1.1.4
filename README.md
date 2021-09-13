# shiny-electron-ngsAnalysis
Electron-Shiny standalone app to bring DNA sequencing data analysis to the point of your mouse-clicks!
## prerequisites
macOS 10.13 (High Sierra) and higher<br/>
R 4.1.1 _(https://cran.r-project.org/bin/macosx/)_<br/>
Python 3.7.7 _(https://www.python.org/downloads/release/python-377/)_<br/>
XQuartz latest _(https://www.xquartz.org/)_<br/>

_you may use the links above to download appropriate versions; double click the .pkg file and follow the on screen instructions to install. XQuartz .pkg is inside the .dmg file._

## install pyhton virtual environment
### open terminal window
> click top menu's Go button and choose utilities<br/>
> inside utilities double-click Terminal<br/>
> it will open up a Terminal window
### check python version
> paste this line below in the terminal prompt and press enter<br/>
`python --version`<br/>

> this will give you the following result<br/>
`Python 2.7.16`<br/>

> this is the mac's system python<br/>
> but we will use python3 and that's why we installed Python 3.7.7, now we will check if the installation was fine<br/>

> copy and paste the following command in your terminal prompt and press enter<br/>
`python3 --version`<br/>

> if everything okay, you will have the following result<br/>
`Python 3.7.7`<br/>

### install virtual environment

> Now paste the following command in your terminal prompt and press enter<br/>
`pip3 install virtualenv`<br/>

> this will install virtual environment and you can confirm the installation by pasting the following command<br/>
`virtualenv`<br/>

> this will produce the following<br/>
`usage: virtualenv [--version] [--with-traceback] [-v | -q] ...` <br/> 

> and the system will exit with error message<br/>
`virtualenv: error: the following arguments are required: dest`<br/>

================================================================<br/>
## install xz (for macOS 10.13)

> The xz library in macOS 10.13 is not compatible with samtools-1.13, which is the ngsAnalyser version.<br/>
>  So, if your macOS is 10.13, you need to install newer xz library.<br/>
>  
>  If you don't have Homebrew installed in your mac already, copy and paste the following command in the terminal prompt, same like you did with virtualenv - Homebrew installation may take few minutes<br/>
>  
>  (following is a single command from ruby upto null, don't miss)<br/>
`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" 2> /dev/null`<br/>

> now use Homebrew to install xz, copy and paste the following command in terminal and press enter<br/>
`brew install xz`<br/>

===============================================================<br/>
## install ngsAnalyser
Download the following file from here<br/> _(https://ifu.cib.csic.es/index.php/s/T7tZ88CwCbY5ECB)_ <br/>
`ngsAnalyser-1.1.4.dmg` <br/>
double click the file and copy the app file to your Applications folder.<br/>

_now you can simply double-click the app icon and it should open a new window on your desktop._

### in case you are missing Command Line Tools
> After successfully completeing all the previous steps, you double-click the ngsAnayser app<br/>

> ngsAnalyser window opens and then gets blurry and nothing works.<br/>

> You can confirm it by clicking the refresh or quit button, there is no response.<br/>

> This is becuase your mac is missing Command Line Tools!<br/>

> The solution is to install Xcode and compile gcc from there.<br/>

> However, if you run the ngsAnalyser with the shell executable, it will run in the terminal.<br/>

> This will directly offer you the prompt to install gcc.<br/>

> Just right click on ngsAnalyser, navigate to _Show Package Contents_ <br/>

> and inside _Contents_ folder, there is a folder named _MacOS_<br/>

> inside _MacOS_, you have the ngsAnalyser shell executable.<br/>

> Double-click it and it will open the app window as well as a terminal window<br/>

> where you have the application log. In absence of gcc, the app fails when trying to install _minipython_<br/>

> Now you will see a auto-prompt asking gcc is missing.<br/>

> Heat the install button and just wait for the installation to finish.<br/>

> Relaunch ngsAnalyser, it should be ready at your service!

## Persistent issues

If you face any further obstackle in running ngsAnalyser, <br/>
please start an issue by going to the issues section above.


### Thanks for using ngsAnalyser!

## copyright

DNA Replication and Genome Integrity Lab, CIB-CSIC<br/>
_(https://www.cib.csic.es/research/cellular-and-molecular-biology/dna-replication-and-genome-integrity)_








