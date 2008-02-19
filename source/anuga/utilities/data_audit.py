"""Track IP of data files in an entire directory tree.
See docstring for the public function IP_verified()
for details.
"""

from os import remove, walk, sep
from os.path import join, splitext

from anuga.utilities.xml_tools import xml2object, XML_element
from anuga.utilities.system_tools import compute_checksum


# Audit exceptions
class NotPublishable(Exception): pass
class FilenameMismatch(Exception): pass
class CRCMismatch(Exception): pass
class Invalid(Exception): pass
class WrongTags(Exception): pass
class Empty(Exception): pass

audit_exceptions = (NotPublishable,
                    FilenameMismatch,
                    CRCMismatch,
                    Invalid,
                    WrongTags,
                    Empty)


def IP_verified(directory,
                extensions_to_ignore=None,
                directories_to_ignore=None,
                files_to_ignore=None,
                verbose=False):
    """Find and audit potential data files that might violate IP

    This is the public function to be used to ascertain that
    all data in the specified directory tree has been audited according
    to the GA data IP tracking process.

    if IP_verified is False:
        # Stop and take remedial action
        ...
    else:
        # Proceed boldly with confidence
        
    verbose controls standard output.
    If verbose is False, only diagnostics about failed audits will appear.
    All files that check OK will pass silently.

    Optional arguments extensions_to_ignore, directories_to_ignore, and
    files_to_ignore are lists of things to skip.

    Examples are:
    extensions_to_ignore = ['.py','.c','.h', '.f'] # Ignore source code
    files_to_ignore = ['README.txt']
    directories_to_ignore = ['.svn', 'misc']

    None is also OK for these parameters.
    
    """

    # Identify data files
    oldpath = None
    all_files = 0
    ok_files = 0
    files_found_in_dir = True
    all_files_accounted_for = True
    for dirpath, filename in identify_datafiles(directory,
                                                extensions_to_ignore,
                                                directories_to_ignore,
                                                files_to_ignore):


        if oldpath != dirpath:
            dir_change = True
            oldpath = dirpath
            files_found_in_dir = False # Reset for this dir
        else:
            dir_change = False

        all_files += 1
        
        basename, ext = splitext(filename)
        license_filename = join(dirpath, basename + '.lic')


        # Look for a XML license file with the .lic
        status = 'OK'
        try:
            fid = open(license_filename)
        except IOError:
            status = 'NO LICENSE FILE'
            all_files_accounted_for = False
        else:
            fid.close()
            
            try:
                license_file_is_valid(license_filename,
                                      filename,
                                      dirpath,
                                      verbose=False)
            except audit_exceptions, e:
                all_files_accounted_for = False                                
                status = 'LICENSE FILE NOT VALID\n'
                status += 'REASON: %s\n' %e

                try:
                    doc = xml2object(license_filename)
                except:
                    status += 'XML file %s could not be read:'\
                              %license_filename
                    fid = open(license_filename)
                    status += fid.read()
                    fid.close()
                else:
                    pass
                    #if verbose is True:
                    #    status += str(doc)


        # Decide if dir header needs to be printed            
        if status != 'OK':
            files_found_in_dir = True
            
                    
        # Only print status if there is a problem (no news is good news)
        if dir_change is True and files_found_in_dir is True:
            print
            print '------------------------------------'
            msg = 'Files without licensing info in dir:'
            print msg, dirpath
            print '------------------------------------'
                

        if status == 'OK':
            ok_files += 1
        else:
            #print dir_change, dirpath, filename + ' (Checksum=%s): '\
            print filename + ' (Checksum=%s): '\
                  %str(compute_checksum(join(dirpath, filename))),\
                  status


    if verbose is True:
        print
        print '---------------------'        
        print 'Audit result for dir: %s:' %directory
        print '---------------------'                
        print 'Number of files audited:  %d' %(all_files)
        print 'Number of files verified: %d' %(ok_files)        
        print

    # Return result        
    return all_files_accounted_for



#------------------
# Private functions
#------------------
def identify_datafiles(root,
                       extensions_to_ignore=None,
                       directories_to_ignore=None,
                       files_to_ignore=None):
    """ Identify files that might contain data

    See function IP_verified() for details about optinoal parmeters
    """

    for dirpath, dirnames, filenames in walk(root):

        for ignore in directories_to_ignore:
            if ignore in dirnames:
                dirnames.remove(ignore)  # don't visit ignored directories


        for filename in filenames:


            # Ignore extensions that need no IP check
            ignore = False
            for ext in extensions_to_ignore:
                if filename.endswith(ext):
                    ignore = True

            if filename in files_to_ignore:
                ignore = True

            if ignore is False:
                yield dirpath, filename


def license_file_is_valid(license_filename, data_filename,
                          dirpath='.', verbose=False):
    """Check that XML license file for given filename_to_verify is valid.

    Input:
        license_filename: XML license file (must be an absolute path name)
        data_filename: The data filename that is being audited
        dir_path: Where the files live
        verbose: Optional verbosity
        

    Check for each datafile listed that

    * Datafile tags are there and match the one specified
    * Fields are non empty (except IP_info which can be left blank)
    * Datafile exists
    * Checksum is correct
    * Datafile is flagged as publishable

    If anything is violated an appropriate exception is raised.
    If everything is honky dory the function will return True.
    """

    doc = xml2object(license_filename)
    
    # Check that file is valid (e.g. all elements there)
    if not doc.has_key('ga_license_file'):
        msg = 'License file %s must have two elements' %license_filename
        msg += ' at the root level. They are\n'
        msg += '  <?xml version="1.0" encoding="iso-8859-1"?>\n'
        msg += '  <ga_license_file>\n'
        msg += 'The second element was found to be %s' %doc.keys()
        raise WrongTags, msg
    

    # Validate elements: metadata, datafile, datafile, ...
    # FIXME (Ole): I'd like this to verified by the parser
    # using a proper DTD template one day....
    # For not, let's check the main ones.
    elements = doc['ga_license_file']
    if not elements.has_key('metadata'):
        msg = 'Tag %s must have the element "metadata"'\
              %doc.keys()[0]
        msg += 'The element found was %s' %elements[0].nodeName
        raise WrongTags, msg

    if not elements.has_key('datafile'):
        msg = 'Tag %s must have the element "datafile"'\
              %doc.keys()[0]
        msg += 'The element found was %s' %elements[0].nodeName
        raise WrongTags, msg    

    for key in elements.keys():
        msg = 'Invalid tag: %s' %key
        if not key in ['metadata', 'datafile']:
            raise WrongTags, msg                    

    
    # Extract information for metadata section
    if verbose: print
    metadata = elements['metadata']

    author = metadata['author']
    if verbose: print 'Author:   ', author
    if author == '':
        msg = 'Missing author'
        raise Exception, msg                
    
    #svn_keywords = metadata['svn_keywords']
    #if verbose: print 'SVN keywords:   ', svn_keywords
    
        
    # Extract information for datafile sections
    datafile = elements['datafile']
    if isinstance(datafile, XML_element):
        datafile = [datafile]


    # Check that filename to verify is listed in license file
    found = False
    for data in datafile:    
        if data['filename'] == data_filename:
            found = True
            break
            
    if not found:
        msg = 'Specified filename to verify %s ' %data_filename
        msg += 'did not appear in license file %s' %license_filename
        raise FilenameMismatch, msg                
            
        
    # Check contents for selected data_filename
    #for data in datafile:
    #    if verbose: print

    # Filename
    if data['filename'] == '':
        msg = 'Missing filename'
        raise FilenameMismatch, msg            
    else:
        filename = join(dirpath, data['filename'])
        if verbose: print 'Filename: "%s"' %filename
        try:
            fid = open(filename, 'r')
        except:
            msg = 'Specified filename %s could not be opened'\
                  %filename
            raise FilenameMismatch, msg

    # CRC
    reported_crc = data['checksum']
    if verbose: print 'Checksum: "%s"' %reported_crc
    
    file_crc = str(compute_checksum(filename))
    if reported_crc != file_crc:
        msg = 'Bad checksum (CRC).\n'
        msg += '  The CRC reported in license file "%s" is "%s"\n'\
               %(license_filename, reported_crc)
        msg += '  The CRC computed from file "%s" is "%s"'\
               %(filename, file_crc)
        raise CRCMismatch, msg
            
    # Accountable
    accountable = data['accountable']
    if verbose: print 'Accountable: "%s"' %accountable
    if accountable == '':
        msg = 'No accountable person specified'
        raise Empty, msg

    # Source
    source = data['source']
    if verbose: print 'Source: "%s"' %source
    if source == '':
        msg = 'No source specified'
        raise Empty, msg                

    # IP owner
    ip_owner = data['IP_owner']
    if verbose: print 'IP owner: "%s"' %ip_owner
    if ip_owner == '':
        msg = 'No IP owner specified'
        raise Empty, msg                                
            
    # IP info
    ip_info = data['IP_info']
    if verbose: print 'IP info: "%s"' %ip_info
    #if ip_info == '':
    #    msg = 'No IP info specified'
    #    raise Empty, msg                                               

    # Publishable
    publishable = data['publishable']
    if verbose: print 'Publishable: "%s"' %publishable
    if publishable == '':
        msg = 'No publishable value specified'
        raise NotPublishable, msg
    
    if publishable.upper() != 'YES':
        msg = 'Data file %s is not flagged as publishable'\
              %fid.name
        raise NotPublishable, msg



    # If we get this far, the license file is OK
    return True
