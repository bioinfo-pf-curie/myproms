#!/bin/bash

### File management
chmod +x ${APPLI_DIR}/cgi-bin/*
mkdir -p ${DATA_DIR}/{exploratory_data,gels,go/goa,go/obo,go/results,logs,peptide_data,quantification_data,spectrum_data,swath_lib,tmp,validation,pathway_data}
if [ ! -e ${APPLI_DIR}/html/data ]; then
    ln -s ${DATA_DIR} ${APPLI_DIR}/html/data
fi
chown -Rh www-data:root ${APPLI_DIR}
chown -Rh www-data:root ${DATA_DIR}

### Update index.html file with ENV variable
sed -i "s@EMAIL_CONTACT@${EMAIL_CONTACT}@g" ${APPLI_DIR}/html/index.html

### Update promsConfig.pm file with ENV variables
if [ ! -e ${APPLI_DIR}/cgi-bin/promsConfig.bck ]; then
    cp ${APPLI_DIR}/cgi-bin/promsConfig.pm ${APPLI_DIR}/cgi-bin/promsConfig.bck
fi
sed -i "s@DB_HOST@${DB_HOST}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@DB_PORT@${DB_PORT}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@DB_USER@${DB_USER}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@DB_PASSWORD@${DB_PASSWORD}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@DB_NAME@${DB_NAME}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@APPLI_DIR@${APPLI_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@DATA_DIR@${DATA_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@SHARED_DIR@${SHARED_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@JAVA_DIR@${JAVA_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@R_DIR@${R_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@QVALITY_DIR@${QVALITY_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@MASSCHROQ_DIR@${MASSCHROQ_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@TPP_DIR@${TPP_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@PYTHON_DIR@${PYTHON_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@MSPROTEO_DIR@${MSPROTEO_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@OPENMS_DIR@${OPENMS_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@PYPROPHET_DIR@${PYPROPHET_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
sed -i "s@HTTP_PROXY@${HTTP_PROXY}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm

### Mascot server
if [ -n "${MASCOT_SERVER_NAME}" ]; then
    sed -i "s@#'MASCOT_SERVER_NAME'@'${MASCOT_SERVER_NAME}'@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
    sed -i "s@MASCOT_SERVER_URL@${MASCOT_SERVER_URL}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
    sed -i "s@MASCOT_SERVER_PROXY@${MASCOT_SERVER_PROXY}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
    sed -i "s@MASCOT_DIR@${MASCOT_DIR}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
    sed -i "s@MASCOT_LINK_FILES_FLAG@${MASCOT_LINK_FILES_FLAG}@g" ${APPLI_DIR}/cgi-bin/promsConfig.pm
fi

### Apache configuration
cat << EOF > /etc/apache2/sites-enabled/000-default.conf
<VirtualHost *:80>

        ErrorLog ${DATA_DIR}/logs/error.log
        CustomLog ${DATA_DIR}/logs/access.log combined

        RewriteEngine On
        RewriteOptions inherit

        DocumentRoot ${APPLI_DIR}/html
        <Directory "${APPLI_DIR}/html">
                Options FollowSymLinks SymLinksIfOwnerMatch
                AllowOverride All
                order allow,deny
                allow from all
                DirectoryIndex index.htm index.html index.php
        </Directory>

        ScriptAlias /cgi-bin/ ${APPLI_DIR}/cgi-bin/
        <Directory "${APPLI_DIR}/cgi-bin">
                Options ExecCGI FollowSymLinks SymLinksIfOwnerMatch
                AllowOverride None
                order allow,deny
                allow from all
        </Directory>

        <IfModule env_module>
                SetEnv PERL5LIB "${APPLI_DIR}/cgi-bin"
        </IfModule>

</VirtualHost>
EOF

cat << EOF >> /etc/apache2/apache2.conf

<Directory ${APPLI_DIR}/>
        Options Indexes FollowSymLinks
        AllowOverride None
        Require all granted
</Directory>
EOF

### Launch Apache
source /etc/apache2/envvars
/usr/sbin/apache2 -DFOREGROUND
