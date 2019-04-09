#!/bin/csh
starver SL18f

set SLEEP = /bin/sleep
set pwddir = $PWD

mkdir -p ${pwddir}/submit

cp submit.xml  ${pwddir}/submit/

cd  ${pwddir}/submit/

star-submit-template -template submit.xml -entities path=${pwddir}

