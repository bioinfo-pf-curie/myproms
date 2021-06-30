# How to connect myProMS to a Mascot server
If you use Mascot search engine, you can configure myProMS to use Mascot search result (.dat) and databanks files directly from your Mascot web server. It is even possible to avoid file duplication between myProMS and Mascot if they share the same file system.

To configure the connection between myProMS and Mascot:
1. **Stop** myProMS server if already running: (`docker-compose down`).
2. If myProMS and Mascot share the **same file system**, edit the file `path/to/myproms/docker-compose.yml` to uncomment lines 49 and 50 and define the pathes to **your** Mascot `data` and `sequence` directories (Only the **left part** of the 2 volume declarations must be edited):
```markdown
   - "/path/to/mascot/data:/mascot/data:ro"
   - "/path/to/mascot/sequence:/mascot/sequence:ro"
```
3. Edit the env file [`path/to/myproms/Mascot/mascot.env`](Mascot/mascot.env) and define the 3 empty variables:
```markdown
MASCOT_SERVER_NAME=<Any name>
MASCOT_SERVER_URL=<http(s)://URL/to/your/Mascot/web/server>
MASCOT_SERVER_PROXY=<Leave empty if no proxy>
```
***Note:*** Use of `https` protocol can sometimes raise security issues. If it is the case, try to use `http` instead.

If myProMS and Mascot are not on the same file system, set the next 2 variables `MASCOT_DATA_DIR` and `MASCOT_SEQUENCE_DIR` to nothing:
```markdown
MASCOT_DATA_DIR=
MASCOT_SEQUENCE_DIR=
```
On the other hand, if they do share the same file system, you can set `MASCOT_LINK_FILES_FLAG=1` if you do not wish to duplicate search result files after import in myProMS. Search result files will be linked instead of copied. 

4. **Copy** the 2 scripts (`myproms4databanks.pl` and `myproms4datFiles.pl`) from the `path/to/myproms/Mascot/` directory to your Mascot `path/to/mascot/cgi` directory. **Make sure** these files are readable and executable by the Mascot server. **Make sure** the Perl shebang of these scripts, `#!/path/to/perl` at line 1 (`#!/usr/local/bin/perl` by default) matches the one of your Mascot installation. You can find it in the  `login.pl` script for instance.

***Security issue:*** You can restrict http access of these 2 scripts to myProMS only by adding the IP address of your myProMS server to the `@remoteHostIPs` Perl array declaration at the beginning of each script (see comments in the scripts). 

5. **Restart** myProMS server: `docker-compose up -d`.


