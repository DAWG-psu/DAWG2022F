# DAWG2022F

### 0. Access to the Roar and open Rstudio server (Be ready to attend the workshop)

#### 0-1. Request Roar Account
#### 0-2. Use Roar - open Rstudio Server
#### 0-3. Install some packages required for the workshop



#### 0-1. Request Roar Account

1. Go to https://www.icds.psu.edu/computing-services/roar-user-guide/
2. Click "Account Request" on top right side of the browser
3. Click "Sign up for your account
4. Log on with your PSU account
5. Fill out the form
    1) Affiliation
    2) University Role
    3) Sponsor Account (Please put your PI's information)
    4) ICDS Linux Clusters (Roar Collab is a new platform that eventually transfer Roar resources in the future, you can mark both of them, but we will use "Roar" still)
  
6. Click "Submit"
7. Wait 1-2 business day and you will get an email to be approved

#### 0-2. Use Roar - open Rstudio Server

1. Go to https://portal2.aci.ics.psu.edu
2. Log on with your PSU account (You need a approved Roar account to access to the portal)
3. Click "Interactive Apps" on the top menu bar

<img width="720" alt="Screen Shot 2022-09-16 at 1 01 09 PM" src="https://user-images.githubusercontent.com/77017866/192337338-fce9edfe-57ee-4d1e-acf0-12c1457bb90c.png">

4. Click "Rstudio Server" on the left menu bar

<img width="720" alt="Screen Shot 2022-09-16 at 1 01 23 PM" src="https://user-images.githubusercontent.com/77017866/192337492-d0a7b951-8e31-4242-ae9a-ad6e832c1d2f.png">

5. Requesting Rstudio server
    1) Allocation - open (unless you already have your own paid allocation from your lab)
    2) DO NOT check "Use custom environment"
    3) DO NOT change anything from "Environment Setup"
    4) Select Version (You would have different method for installing packages based on the version you select - you can select the latest one, but DO NOT select v3.x.x, it might be too old for some packages
    5) Select Number of hours (up to 48 hours - if you request more times, it will take more time to open the session, so do not request more than what you need)
    6) Select Number of Cores (up to 20 - same has number of hours, amouont of request resources will affect the time to start the session)
    7) Select Memory per core (in GB) (TOTAL memory you can request in 128 gb ( Number of cores * memory per core = total memory)
    8) Select Note type (If you are using "Open" allocation, just leave this as "ACI-b Standard Core")
    9) Click "Launch" at the bottom
    
6. Wait until the session is open - it would take serveral minutes, if your session does not open later 30 minutes, there might be something wrong with your setting. 

<When you just launch your session and still "Queued">
<img width="720" alt="Screen Shot 2022-09-16 at 1 01 58 PM" src="https://user-images.githubusercontent.com/77017866/192339620-4d97f394-aab4-430d-9254-92dd312bb821.png">

<After several minutes, you session is ready to run ("Running")>
<img width="720" alt="Screen Shot 2022-09-16 at 1 12 49 PM" src="https://user-images.githubusercontent.com/77017866/192339650-ad06199e-3257-4965-951b-52375014be1a.png">

7. Click "Connect to Rstudio Server
    
#### 0-3. Install some packages required for the workshop

1. Install BioConductor (https://bioconductor.org/install/)
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
```
* In the last code version = "3.14" should be changed based on your R version you requested from the previous section based on the table attached here

<img width="273" alt="Screen Shot 2022-09-16 at 1 13 10 PM" src="https://user-images.githubusercontent.com/77017866/192342143-ce72035c-cd85-477c-8312-af9c89655a74.png">

For example, if you open R (4.1.x) you need to install Bioconductor v 3.13 or 3.14, if you open R (4.0.x), you need to install Bioconductor v3.12. or v3.11.

IF you are using your own machine, and your R version is higher than v 4.2.0, you need to install Bioconductor v3.15

2. Install DADA2 and their dependencies (https://benjjneb.github.io/dada2/dada-installation.html)
(From now on, I am following with R v4.1.2 and Bioconductor v3.14)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.14") 
```

* If you are getting error message saying "package(s) not installed when version(s) same as current; use `force = TRUE` to re-install:", you already have installed either or both of the packages, you can check your installation by typing
```
library("dada2")
library("ShortRead")
```

Now you are ready to come to the workshop, before the workshop, I will upload files and script here and let you know. 


