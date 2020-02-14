# myProMS Server
Mass spectrometry(MS)-based proteomics data management and analysis software.
**myProMS** is a web server developed to allow MS specialists and their collaborators to efficiently manage, process and interprete MS-based proteomics data.
Go to **myProMS** [home page](http://myproms-demo.curie.fr) for more information about the software and access an online **demo** version.
myProMS can run on any operating system compatible with [Docker](https://www.docker.com).

## Version
Current version is 3.9.1.
Check [releases](https://github.com/bioinfo-pf-curie/myproms/releases) for versions history.

## License
myProMS is freely available under [CeCILL license](LICENSE).

## Requirements
* Docker: [https://www.docker.com](https://www.docker.com)
Two Docker images are required:
  * [myproms/myproms\_base:1.3.2](https://hub.docker.com/r/myproms/myproms_base): Tools and system dependencies.
  * [mysql:5.5](https://hub.docker.com/_/mysql): For database management).
These images will be automatically pulled at first launch of the server.
* Git: [https://git-scm.com](https://git-scm.com) - *Recommanded*

## Repository
This repository contains:
* All HTLM/JavaScript/Image/MySQL/Perl/R/Java/bash files/scripts specific to myProMS Server.
* Scripts, configuration files and instructions for installing and running the application.

## Installation and run instructions
**Note**: You might want to read the [Custom configuration](#custom-configuration) section before to decide whether the default configuration suits you.
1. Make sure docker is running.
2. Open a console.
3. Pull the required docker images (optional):
```
docker pull myproms/myproms_base:1.3.2
docker pull mysql:5.5
```
4. Create a home directory for your docker projects (optional):
```
mkdir docker
cd docker
```
5. Retrieve myproms **repository** from GitHub:  
`git clone --depth 1 https://github.com/bioinfo-pf-curie/myproms myproms` - *Recommanded*  
Or download and extract a compressed archive of the repository (https://github.com/bioinfo-pf-curie/myproms).  
This will create a directory named `myproms` containing all files on this repository.  
**Important**:
    * **Linux/MacOS** users: Make sure the file `path/to/myproms/start.sh` is executable (`ls -l path/to/myproms/start.sh`). If not, run the command `chmod +x path/to/myproms/start.sh`.
    * **Windows** users: Make sure the file `path/to/myproms/start.sh` has Unix-style line endings `(\n)`. It should be the case by default unless the file has been copied/altered.  
6. To **start** myProMS Server with **docker-compose**:
```bash
cd to/directory/myproms
docker-compose up -d
```
Two containers named `myproms-mysql` and `myproms-server` should be running. You can check this with the command: `docker ps`.

7. To **use** myProMS:  
Open a web browser and go to `http://localhost:8080` or `http://127.0.0.1:8080` to access the server home page. myProMS default **user account** is `login` with **password** `password`. Once connected, go to the **Users** section to create at least one *real* `bioinformatician` account. Logout and re-connect using this new account.
We strongly recommand you then close the `login` account or at least change its password.  
8. To **shut down** myProMS Server:  
From `myproms` directory type `docker-compose down`  
9. Repeat steps **6 to 8** to use the application as often as you wish.  

## Custom configuration
Some configuration parameters can be customized to match your preferences. In particular, we recommand to at least change the default MySQL user password.

#### Database sensitive parameters
Sensitive database variables are strored in the files [myproms-mysql.env](myproms-mysql.env) and [myproms-server.env](myproms-server.env). Open these files with any text editor and make changes as indicated below.  
The database **user** and **password** are defined in the following paired variables:  
* **User**: `MYSQL_USER` and `DB_USER` (default: `myproms`)  
* **password**: `MYSQL_PASSWORD` and `DB_PASSWORD` (default: `myproms`)

**Important**: If you have started myProMS Server at least once before making any configuration changes, you must reset myProMS configuration file to it's original state for your changes to be applied at the next start. This is done with the following command:  
`cp path/to/myproms/myproms_appli/cgi-bin/promsConfig.bck path/to/myproms/myproms_appli/cgi-bin/promsConfig.pm`  
The database connection credentials will be transferred to the `promsConfig.pm` file. Make sure user access to this file and the two `.env` files is under tight control.

#### Time zone
myProMS time zone is set to `Europe/Paris` by default. This can be changed by editing the variable `TZ` in the two `.env` files described above.

#### Network ports
The default port used to access the myProMS server with your browser is set to `8080`. You can change this by opening the file [docker-compose.yml](docker-compose.yml) with any text editor and go to the section below:
```markdown
 ports:
   - "8080:80"
```
**Note**: You can use the default Apache port `80` (`   - "80:80"`) if no other services installed on your system are served on it. This will allow you to skip the port declaration in myProMS URL: `http://localhost`.  

#### Your email contact
An email contact is displayed on myProMS startup page. The default is `myproms@curie.fr`. You recommand using a contact local to your institute instead.

## Data persistence
Because Docker containers are created/deleted each time you start/shut down the myProMS Server, all data are stored outside these containers through mounted [volumes](https://docs.docker.com/storage/volumes/) to insure data persistence.
* **MySQL database** data are stored in the directory `path/to/myproms/myproms_db`
* **Flat file** data are stored in the directory `path/to/myproms/myproms_data`

## Contact
**Email**: myproms@curie.fr

