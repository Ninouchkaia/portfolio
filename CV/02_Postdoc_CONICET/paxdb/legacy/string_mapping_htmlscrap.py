#!/usr/bin/python
# -*- coding: utf-8 -*-



import timeit
start = timeit.default_timer()
import sys
import os
from commands import getoutput # permet d'obtenir l'output d'une commande bash.



#liste = open("/home/nina/text_files/%s" % sys.argv[1])
##for example : 'no_uniprot_match.txt'
#string_ids = liste.readlines()
#liste.close()


#for string_id in string_ids:
#string_id = string_id.replace('\n', '') # remove '\n' only
#print string_id
#cmd1 = "wget -w 5 \"http://pax-db.org/#!protein/1110761\" -O 1110761.html " 
cmd2 =  "curl \"http://pax-db.org/#!search?q=10090.ENSMUSP00000000940\" -L -o dumpfile -w 'Last URL was: %{url_effective}'"
#% (string_id, string_id) # pour ubuntu
os.system(cmd2)
#liste.close()



stop = timeit.default_timer()
print stop - start 