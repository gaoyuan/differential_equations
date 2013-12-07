Hangman Solver
==========
An AI to play the game of Hangman and save humans.

Environment
------------
This application was built and tested on Mac OS X.

Dependencies
------------
1. [Node.js](http://nodejs.org/)

2. The request module

	sudo npm install -g request

3. [Coffeescript](http://coffeescript.org/)

	sudo npm install -g coffee-script

Usage
------------
$coffee main.coffee

To visualize the logs, start a local http server under the log directory.
For example if you have python installed,

$python -m SimpleHTTPServer 8080

And type in your favrorite browser,

http://localhost:8080/visualize.html

Acknowledgement
------------
The word frequency data used by the AI is from http://crr.ugent.be/papers/SUBTLEX-UK.txt
