"""Validate_all.py 




"""

import os, time, sys


# List any sub directory to exclude from validation.
# Current working directory ('.') should always be excluded to avoid 
#infinite recursion
dirs_to_skip = ['.'] # Always skip current dir
#dirs_to_skip += ['patong_beach_validation'] # This takes about 40h

validation_dirs_and_files = []
for dirpath, dirnames, filenames in os.walk('.'):

    if '.svn' in dirnames:
        dirnames.remove('.svn')  # don't visit SVN directories

    dir = os.path.split(dirpath)[-1]
    if dir in dirs_to_skip:
        #print 'Skipping %s' % dirpath
        continue
    
    #print 'Searching dir', dirpath
   

    for filename in filenames:
        if filename.startswith('validate_') and filename.endswith('.py'):
            #print 'Found %s in %s' %(filename, dirpath)
            validation_dirs_and_files.append((dirpath, filename))            

    # get repeatable order on different machines
    validation_dirs_and_files.sort()


print 
print '---------------------------------------------------------'
print 'Running all validation tests - some may take many minutes'
print 'and some may require memory in the order of 8-16GB       '
print '---------------------------------------------------------'

print 'Validation test suites:'
for path, filename in validation_dirs_and_files:
     print '    ', os.path.join(path, filename)
print
print

t0 = time.time()
parentdir = os.getcwd()


for path, filename in validation_dirs_and_files:
    # print 'filename path', path, filename

    os.chdir(path)
    s = 'python %s' % filename
    print 60*'-'
    print s
    print 60*'-'
    os.system(s)
    
    # Back to parent directory
    os.chdir(parentdir)

    # print 'current dir', os.getcwd() 
    
print 'That took %.2f seconds in total' %(time.time()-t0)
    

        