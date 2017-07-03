CSB=`python -c "import csb, os; print(os.path.dirname(csb.__path__[0]))"`

epydoc --html -v -o /tmp --name CSB --no-private --introspect-only --exclude csb.test.cases --css $CSB/epydoc.css --fail-on-error --fail-on-warning --fail-on-docstring-warning csb
python $CSB/csb/test/app.py

read -p "Press ENTER to exit..." NULL
