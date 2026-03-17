# TISC Frequently Asked Questions

### Q1. How can I learn to use this code on my own?
**A1.** After download, follow instructions in `tisc/doc/README` file in order to compile. Then go to `tisc/demo` and try the examples there to ensure that the code is working properly. Then read carefully the documentation and `tisc/doc/template.PRM` and try to experiment modifying the `*.PRM` files of the demos. Then learn how to build `*.UNIT` input files from the examples and using the explanations. Make sure you did all this before contacting the author!

### Q2. I have compilation errors, is there a problem in the program?
**A2.** I tested TISC only in some UNIX operative systems such as Linux (Ubuntu), Mac OSX 10, and Solaris (from SUN). Most probably you have a problem with the compilation options or with the libraries. Try to ask help among your colleagues to succeed in the compilation.

### Q3. I compiled the program, but I never managed to run it, and I always get 'segmentation fault' errors.
**A3.** This probably means you have some incompatibility between some of the functions called in TISC and the ones available in the libraries of your operative system. You have to edit the code and search for those problems. You will need the help of someone around you who has experience in programing or compiling C programs.

### Q4. Can you help me to solve the problems I have to run TISC?
**A4.** In general, I try to answer the questions you post on GitHUb, but often i cannot say much more than what you just read in this file. Technical support on the compilation process is not provided. Because my interests are in the scientific use of TISC and its application to geological scenarios, I will be pleased to help you with any problem you have during the development of a model and will be pleased to establish a scientific collaboration.

### Q5. Is it possible to stop a fault and then reactivate it again?
**A5.** You can stop a fault at a given time by using the `time_stop` parameter in the `*.UNIT` file. To reactivate it, repeat the same fault geometry in another `*.UNIT` file with the velocity you wish. See `doc/template.UNIT` to get the possible parameters in every `*.UNIT` file.