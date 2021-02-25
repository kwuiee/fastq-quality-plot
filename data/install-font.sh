set -e

tar zxvf wqy-zenhei.tar.gz
cp -rf wqy-zenhei /usr/share/fonts/wqy-zenhei
cp -f 65-wqy-zenhei.conf /etc/fonts/conf.d/
rm -rf wqy-zenhei 65-wqy-zenhei.conf
