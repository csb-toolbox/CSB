# USAGE: build.sh [pythonx.y]
#    pythonx.y  Python interpreter executable (default="python")

# In cron this script must be executed in this way:
#   bash -lc 'build.sh'
# otherwise env vars from ~/.bashrc profile may not be initialized ($PYTHONPATH)

PLATFORM=''
PYTHON='python2'
if [ "$1" ]; then
	PYTHON=$1
	PLATFORM=`$PYTHON -c "import sys; print(sys.version_info.major)"`
fi

ROOT=$HOME/BUILD									# Build directory - please configure
CHECKOUT=$ROOT/$PYTHON
OUTPUT=$CHECKOUT/CSB/
LOG=$OUTPUT/log

SMTP='localhost'									# SMTP server (e.g. smtp.gmail.com for Gmail or localhost for a local SMTP)
BOTMAIL='csb.build@domain.com'						# this will be the FROM: email field - please configure
BOTPWD=''											# password for $BOTMAIL, needed only if $SMTP requires authentication (e.g. Gmail)
OPERATORS=$BOTMAIL									# this will be the TO: email field - please configure
LOGURL="http://`hostname -f`/csb$PLATFORM/log"

mkdir -p $ROOT
mkdir -p $CHECKOUT

# Checkout
cd $CHECKOUT
rm -f -r CSB
git clone https://github.com/csb-toolbox/CSB.git

# Build
cd $OUTPUT 
$PYTHON $OUTPUT/tip/csb/build.py -o $OUTPUT -v 1 >> $LOG 2>&1

# Zip in order to wrap in a static file name: build.zip
zip -r build.zip csb-*.tar.gz 2>&1

# Clean
rm -f -r $OUTPUT/tip

# Prepare email
LOGBODY=`python -c "print open('$LOG').read().replace('\n', '<br>').replace('\'', '').replace('\"', '')"`
VERSION=`cat $LOG | grep "# Done (" | python -c "import sys; print sys.stdin.read()[8:19].split(')')[0]"`
FAILED=`grep -m 1 "  DID NOT PASS" $LOG`
EX=0

if [ "$VERSION" ]
	then
		if [ "$FAILED" ]
			then	
				MAIL="Subject: Broken CSB Build: $VERSION\n\n$VERSION: $FAILED\n\nBuild log:\n\n $LOGBODY\n"
				EX=99
			else
				PASSED=`grep -E "  Passed all .+ tests" $LOG | tail -1`
				MAIL="Subject: Normal CSB Build: $VERSION\n\n$VERSION: $PASSED\n\nBuild log: $LOGURL\n"
		fi
	else
		MAIL="Subject: Broken CSB Build: CRASH\n\nThe build bot has crashed. Build log:\n\n $LOGBODY\n"
		EX=98
fi

# Notify recipients
if [ "$BOTPWD" ]
	then
		$PYTHON -c "import smtplib; s = smtplib.SMTP_SSL('$SMTP'); s.login('$BOTMAIL', '$BOTPWD'); s.sendmail('$BOTMAIL', '$OPERATORS', '$MAIL'.replace('<br>', '\n'));"
	else
		$PYTHON -c "import smtplib; smtplib.SMTP('$SMTP').sendmail('$BOTMAIL', '$OPERATORS', '$MAIL'.replace('<br>', '\n'));"
fi

# Done
exit $EX

