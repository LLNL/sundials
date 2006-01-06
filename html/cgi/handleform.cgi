#!/usr/bin/perl

###########################################################################
#
# Based on:
#
# FormHandler                       Version 2.0
# Written by Matthew Wright         mattw@worldwidemart.com
# Created 05/31/96                  Last Modified 3/6/97
#
# parse_form()                      Version 1.5
# Written by Matthew Wright         mattw@worldwidemart.com
# Created 9/30/96                   Last Modified 3/28/97
#
# email_check()                     Version 1.0
# Written by Matthew Wright         mattw@worldwidemart.com 
# Created 8/1/96                    Last Modified 3/23/97 
#
# send_email()                      Version 1.0
# Written by Craig Patchett         craig@patchett.com
#     and Matthew Wright            mattw@worldwidemart.com
# Created 11/16/97                  Last Modified 11/16/97
#
# parse_template()                  Version 2.0
# Written by Craig Patchett         craig@patchett.com
#    and Matthew Wright             mattw@worldwidemart.com
# Created 9/9/96                    Last Modified 4/6/97
#
# format_date()                     Version 1.5
# Written by Craig Patchett         craig@patchett.com
# Created 8/15/96                   Last Modified 5/31/97
#
# lock()                            Version 2.0
# Written by Craig Patchett         craig@patchett.com
# Created 9/16/96                   Last Modified 3/28/97 
#
# unlock()                          Version 2.0
# Written by Craig Patchett         craig@patchett.com
# Created 9/16/96                   Last Modified 3/28/97
#
# Copyright 1997 Craig Patchett & Matthew Wright.  All Rights Reserved.
# This program is part of The CGI/Perl Cookbook from John Wiley & Sons.
# License to use this program or install it on a server (in original or
# modified form) is granted only to those who have purchased a copy of The
# CGI/Perl Cookbook. (This notice must remain as part of the source code.)
#
###########################################################################
#
# Modified by Radu Serban for SUNDIALS.
# January 14, 2005.
#
###########################################################################

#--------------------------------------------------------------------------
# Define configuration constants
#--------------------------------------------------------------------------

# @REFERERS contains the host names and IP addresses of the domains which
# are allowed to use this copy of FormHandler to parse their forms.

@REFERERS = ('localhost:8090');
#@REFERERS = ('www.llnl.gov');
#@REFERERS = ('cmg.llnl.gov');

# %CONFIG defines which form fields should be considered configuration
# fields rather than standard data fields. Each of the default variables
# defined in the array below have special meaning to FormHandler and are
# usually set using hidden fields. Default values used in the array will
# be overridden by form fields with the same name. Any variable that should
# be considered a configuration variable must be defined in this array.

%CONFIG = ('required',                '',
           'realname',                '', 
           'email',                   '',
           'field_names',             '',
           # Stuff for registration
           'subject',                 '', 
           'redirect',                '',
#           'developer_addr',          'alanh@llnl.gov,cswoodward@llnl.gov,radu@llnl.gov',
           'developer_addr',          'radu@llnl.gov',
           'developer_names',         'Carol, Radu, and Alan',
           'email_dl_template',       'email_developers.txt',
           'log_dl_template',         'log_dl_templ.html', 
           'log_dl_filename',         'logfile_dl.txt',
           # Stuff for subscription
           'action_sundials-announce','',  
           'action_sundials-users',   '',
           'log_action',              '',
           'majordomo_addr',          'majordomo@lists.llnl.gov',
           'subject_ml',              'Mailing List Request', 
           'email_ml_template',       'email_ml.txt',
           'log_ml_template',         'log_ml_templ.html', 
           'log_ml_filename',         'logfile_ml.txt',
           # error and success templates
           'error_template',          'missing.html',
           'success_html_template',   'confirm.html');

# Host name of the web server.

$WEB_SERVER = 'localhost:8090';
#$WEB_SERVER = 'www.llnl.gov';
#$WEB_SERVER = 'cmg.llnl.gov';

# sendmail program

$SENDMAIL = '/usr/sbin/sendmail';

# Directory where all lock files will be placed.

$LOCK_DIR = 'lock/';

# Max number of seconds the lock script will wait before overiding the lock file.

$MAX_WAIT = 5;

# Directory in which all of your required files are placed.

$REQUIRE_DIR = './';

############################################################################
# Get Required subroutines which need to be included.                      #
############################################################################

# Push the $REQUIRE_DIR onto the @INC array for include file directories.

push(@INC, $REQUIRE_DIR) if $REQUIRE_DIR;

# Require Necessary Routines for this script to run.

#--------------------------------------------------------------------------
# Check that the form is coming from a web site that is included in the
# @REFERERS array.  If not, sent out an error. Otherwise, continue.
#--------------------------------------------------------------------------

# Set the flag to 0 so the referer will fail by default.

$check_referer = "0";

# Get the hostname out of the referring URL.  If there is no 
# referring URL, set flag to 1, since some browsers don't pass the 
# HTTP_REFERER environment variable.

if ($ENV{'HTTP_REFERER'} =~ m#http://([^/]+)#) {
    $hostname = $1;
}
else {
    $check_referer = 1;
}

# If there is a hostname found, check it against the hostnames in the 
# @REFERERS array until a match is found.  Set flag to 1 if match is 
# found.

if ($hostname) {
    foreach $referer (@REFERERS) {
        if ($hostname =~ /$referer/i) {
            $check_referer = 1;
            last;
        }
    }
}

# If flag is not set to 1, throw an error to the &error subroutine.

if (!$check_referer) {
    &error("$ENV{'HTTP_REFERER'} is not allowed access to this program.");
}

#--------------------------------------------------------------------------
# Get dates for use in templates and logs.
#--------------------------------------------------------------------------

# Get the current time in seconds.

$current_time = time;

# Create a log_date for use in the logs using normal clock and date format.

$CONFIG{'log_date'} = &format_date($current_time, "<mtimesec> <0m>/<0d>/<yr>");

#--------------------------------------------------------------------------
# Parse the form contents and put configuration fields into %CONFIG and    
# other fields in %FORM.  This is an external subroutine which was required
#--------------------------------------------------------------------------

if (!&parse_form) {
    &error($Error_Message);
}

#--------------------------------------------------------------------------
# If the field names configuration field has been filled in, parse it and
# define field names which may be used later in the script.
#--------------------------------------------------------------------------

if ($CONFIG{'field_names'}) {
    @field_names = split(/&/, $CONFIG{'field_names'});

    # For each of the field names specified in the field_names form 
    # field, split the old name and new name by = sign and then set the new 
    # names in the %ALT_NAME array.
    
    foreach $field_name (@field_names) {
        ($def_name, $def_value) = split(/=/, $field_name);
        $ALT_NAME{$def_name} = $def_value;
    }
}

#--------------------------------------------------------------------------
# Check the email field and if it contains an obviously bad e-mail address.
#--------------------------------------------------------------------------

if (!&email_check($CONFIG{'email'})) {
    $CONFIG{'email'} = "";
}

#--------------------------------------------------------------------------
# Check any fields which were required.
#--------------------------------------------------------------------------

# If this is the mailing list subscription,
# require that at least one of the mailing lists be chosen.

if ($FORM{'subscribe'}) {
    if (!($FORM{'sundials-announce'}) && !($FORM{'sundials-users'})) {
        push(@missing_required_fields, "Your mailing list selection(s).");
    }
}

# Split required fields by commas in the 'required' configuration field.

@required = split(/,/, $CONFIG{'required'});

# For each of the form fields specified in the @required array, 
# determined from the required form field, make sure it is not 
# blank.  Push any blank fields onto the @missing_required_fields array.

foreach $required (@required) {
    if (!($CONFIG{$required}) && !($FORM{$required})) {
        push(@missing_required_fields, $required);
    }
}

# If the @missing_required_fields array exists, then thrown an error.

if (@missing_required_fields) {
    &error('missing_required_fields');
}

#--------------------------------------------------------------------------
# DOWNLOAD REGISTRATION
#--------------------------------------------------------------------------

if ( $FORM{'register'} ) {

    # Log registration info to file

    &log2file($CONFIG{'log_dl_template'}, $CONFIG{'log_dl_filename'});
    
    # Send notification to developers

    $from = "($CONFIG{'realname'}) $CONFIG{'email'}";

    if (&send_email($CONFIG{'subject'}, $from, $CONFIG{'developer_addr'}, 
                    '', '',
                    $CONFIG{'email_dl_template'})) {
        &error($Error_Message);
    }

    # Test if subscription to a mail list was desired

    if ( $FORM{'sundials-announce'}  || $FORM{'sundials-users'} ) {
        $action = 'subscribe';
    } else {
        $action = '';
    }
}

if ( $FORM{'subscribe'} ) {
    $action = $FORM{'action'};
}

#--------------------------------------------------------------------------
# MAILING LIST REGISTRATION
#--------------------------------------------------------------------------

if ( $action ) {

    $space = " ";

    # Set-up log_action for logging

    $CONFIG{"log_action"} = $action;
    if ($FORM{'sundials-announce'}) {
        $CONFIG{"log_action"} .= $space;
        $CONFIG{"log_action"} .= $FORM{'sundials-announce'};
    }
    if ($FORM{'sundials-users'}) {
        $CONFIG{"log_action"} .= $space;
        $CONFIG{"log_action"} .= $FORM{'sundials-users'};
    }

    # Log registration info to file

    &log2file($CONFIG{'log_ml_template'}, $CONFIG{'log_ml_filename'});

    # Set-up variables for email to majordomo

    if ( $FORM{'sundials-announce'} ) {
        $CONFIG{"action_sundials-announce"} = $action;
        $CONFIG{"action_sundials-announce"} .= $space;
        $CONFIG{"action_sundials-announce"} .= $FORM{'sundials-announce'};
        $CONFIG{"action_sundials-announce"} .= $space;
        $CONFIG{"action_sundials-announce"} .= $CONFIG{'email'};
    }

    if ( $FORM{'sundials-users'} ) {
        $CONFIG{"action_sundials-users"} = $action;
        $CONFIG{"action_sundials-users"} .= $space;
        $CONFIG{"action_sundials-users"} .= $FORM{'sundials-users'};
        $CONFIG{"action_sundials-users"} .= $space;
        $CONFIG{"action_sundials-users"} .= $CONFIG{'email'};
    }

    # Send request to majordomo

    if (&send_email($CONFIG{'subject_ml'}, $CONFIG{'email'}, $CONFIG{'majordomo_addr'}, 
                    '',  '', 
                    $CONFIG{'email_ml_template'})) {
        &error($Error_Message);
    }

}

#--------------------------------------------------------------------------
# Send back HTML response.  
#--------------------------------------------------------------------------

if ( $FORM{'register'} ) {

    # Set a cookie valid for 1 year and then redirect.
    
    $exp_time=gmtime(time()+365*24*3600)." GMT";  # Add 12 months (365 days)
    print "Set-cookie: SunRegCookie=OK; path=/; expires=$exp_time;\n";
    print "Location: $CONFIG{'redirect'}\n\n";

}

if ( $FORM{'subscribe'} ) {
    
    # Confirm subscription

    print "Content-type: text/html\n\n";
    if (!&parse_template($CONFIG{'success_html_template'}, *STDOUT)) {
        &error("Can't open $CONFIG{'success_html_template'} ($!).");
    }

}

exit;

###########################################################################
#
#  LOCAL SUBROUTINES
#
###########################################################################

#--------------------------------------------------------------------------
# Log to file, using template
#--------------------------------------------------------------------------

sub log2file {
    
    my $log_template = $_[0];
    my $log_filename = $_[1];
    
    # Lock the log file and open it for appending.

    if (&lock($log_filename, $LOCK_DIR, $MAX_WAIT)) {
        &error($Error_Message);
    }
    open(LOG, ">>$log_filename")
        || &error("Could not open $log_filename ($!).");
    
    # Log data to file using template

    &parse_template($log_template, *LOG)
        || &error("Could not open $log_template ($!).'}");
    
    # Close and unlock the log file
    
    close(LOG);
    &unlock($log_filename, $LOCK_DIR);

}

#--------------------------------------------------------------------------
# Print out the Error Messages and exit.
#--------------------------------------------------------------------------

sub error {

    # Localize the error variable.
    
    local($error) = $_[0];
    print "Content-type: text/html\n\n";
    
    # If the error is because of missing_required_fields, use the
    # error_template that was specified or print out a generic response.
    
    if ($error eq 'missing_required_fields') {
    
        # Prepare the error_fields config field so that users can use it 
        # in their templates.
        
        $CONFIG{'error_fields'} = "<UL>\n";
        foreach $missing_required_field (@missing_required_fields) {
            if ($ALT_NAME{$missing_required_field}) {
                $CONFIG{'error_fields'} .= "<LI>$ALT_NAME{$missing_required_field}</LI>\n";
            }
            else {
                $CONFIG{'error_fields'} .= "<LI>$missing_required_field</li>\n";
            }
        }
        $CONFIG{'error_fields'} .= "</UL>";
        
        # Print out formatted template to user.
        
        if (!&parse_template($CONFIG{'error_template'}, *STDOUT)) {
            $error = "Can't open $CONFIG{'error_template'} ($!).";
        }
        else { exit }
        
    }
    
    # For any other errors, just print a title and heading which supplies 
    # the error.
    
    print <<HTML_END;
<HTML>
   <HEAD>
      <TITLE>$error</TITLE>
   </HEAD>
   <BODY BGCOLOR=#FFFFFF TEXT=#000000>
      <CENTER>
      <H4>$error</H4>
      </CENTER>
   </BODY>
</HTML>
HTML_END
    exit;
}

#--------------------------------------------------------------------------
# parse_form
#--------------------------------------------------------------------------
#
# Function:      Takes form field data from a POST or GET request and      #
#                converts it to name/value pairs in the %FORM array or,    #
#                if a corresponding entry in the %CONFIG array is defined, #
#                in %CONFIG.                                               #
#                                                                          #
# Usage:         &parse_form;                                              #
#                                                                          #
# Variables:     None                                                      #
#                                                                          #
# Returns:       0 if invalid request method                               #
#                1 if successful                                           #
#                                                                          #
# Uses Globals:  Sets %CONFIG with name/value pairs if corresponding entry #
#                  in %CONFIG is defined                                   #
#                Otherwise sets entries in %FORM                           # 
#                $Error_Message for descriptive error messages
#

sub parse_form {
    local($name, $value, $pair, $buffer, @pairs);
    
    # Check for request method and handle appropriately
    
    if ($ENV{'REQUEST_METHOD'} eq 'GET') {
        @pairs = split(/&/, $ENV{'QUERY_STRING'});
    }
    elsif ($ENV{'REQUEST_METHOD'} eq 'POST') {
        read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
        @pairs = split(/&/, $buffer);
    }
    else {
        $Error_Message = "Bad request method ($ENV{'REQUEST_METHOD'}).  Use POST or GET";
       return(0);
    }

    # Convert the data to its original format
    
    foreach $pair (@pairs) {
        ($name, $value) = split(/=/, $pair);

        $name =~ tr/+/ /;
        $name =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
        $name =~ s/\n//g;
        $value =~ tr/+/ /;
        $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
        $value =~ s/\n//g;

        # If they try to include server side includes, erase them, so they
        # arent a security risk if the HTML gets returned.  Another
        # security hole plugged up.
        
        $value =~ s/<!--(.|\n)*-->//g;

        # Store name/value pair in %CONFIG if the corresponding entry is
        # defined
        
        if ($CONFIG{$name}) {
            $CONFIG{$name} .= ",$value";
        }
        elsif (defined($CONFIG{$name})) {
            $CONFIG{$name} = $value;
        }
        
        # Otherwise store in %FORM
        
        elsif ($FORM{$name}) {
            $FORM{$name} .= ",$value";
        }
        else {
            $FORM{$name} = $value;
        }
    }
    return(1);
}

#--------------------------------------------------------------------------
# Simple email syntax check
#--------------------------------------------------------------------------
#
# Function:      Checks an email address to see if it passes a simple      #
#                syntax check. (This routine will not check to see if the  #
#                address is an actual address.)                            #
#                                                                          #
# Usage:         &email_check($email_address);                             #
#                                                                          #
# Variables:     $email_address -- String containing the address to check  #
#                                  Example: 'someone@somewhere.com'        #
#                                                                          #
# Returns:       0 if the email address is invalid                         #
#                1 if the address is in a valid format                     #
#

sub email_check {
    local($email) = $_[0];

    # Check that the email address doesn't have 2 @ signs, a .., a @., a 
    # .@ or begin or end with a .

    # Allow anything before the @, but only letters numbers, dashes and 
    # periods after it.  Also check to make sure the address ends in 2 or 
    # three letters after a period and allow for it to be enclosed in [] 
    # such as [164.104.50.1]

    if ($email =~ /(@.*@)|(\.\.)|(@\.)|(\.@)|(^\.)|(\.$)/ || 
        ($email !~ /^.+\@localhost$/ && 
         $email !~ /^.+\@\[?(\w|[-.])+\.[a-zA-Z]{2,3}|[0-9]{1,3}\]?$/)) {
        return(0);
    }

    # If it passed the above test, it is valid.
    
    else {
        return(1);
    }
}

#--------------------------------------------------------------------------
# Sends an email message via  a pipe connection to sendmail
#--------------------------------------------------------------------------
#                                                                          #
# Function:      Sends an email message and optional attached files via    #
#                a pipe connection to sendmail or a compatible program.    #
#                                                                          #
# Usage:         &send_email($subject, $from, $to[, $cc, $bcc, $body,      #
#                            $files, $encoding]);                          #
#                                                                          #
# Variables:     $subject --  String containing subject of message.        #
#                             Example 'Buy the CGI/Perl Cookbook!'         #
#                $from --     String containing email address of person    #
#                             sending message. An associated name can      #
#                             follow the address if placed in parentheses. #
#                             Example 'me@home.com (My Name)'              #
#                $to --       String containing email addresses to send    #
#                             message to. Multiple addresses should be     #
#                             separated by commas. Associated names        #
#                             can follow each address if placed in         #
#                             parentheses.                                 #
#                             Example 'him@place.com (Name),her@place.com' #
#                $cc --       String containing email addresses to send    #
#                             copies of the message to. Same format as $to.#
#                $bcc --      String containing email addresses to send    #
#                             blind copies of the message to (i.e., nobody #
#                             else receiving the message will know that    #
#                             copies were sent to these addresses). Same   #
#                             format as $to.                               #
#                $body --     Full path to file containing the body of the #
#                             message or text containing body of message   #
#                             (if $body doesn't begin with a directory     #
#                             delimiter and contains at least one space    #
#                             then the subroutine assumes it contains      #
#                             message text).                               #
#                             Example '/home/user/body.txt'                #
#                             Example 'This is message text.'              #
#                                                                          #
# Returns:       Nothing                                                   #
#                                                                          #
# Uses Global:   $SENDMAIL for the path to the sendmail program (assumes   #
#                          that sendmail is in search path if not set)     #
#

sub send_email {

    local($subject, $from, $to, $cc, $bcc, $body) = @_;
    local($i, $name, $status, $message) = '';

    # Attempt to set default values if globals aren't set
    
    if (!$SENDMAIL) { $SENDMAIL = "sendmail" }
    
    # Split the input into arrays where needed, since values are passed
    # as strings separated by commas.

    local(@to) = split(/, */, $to);
    local(@cc) = split(/, */, $cc);
    local(@bcc) = split(/, */, $bcc);

    # Set up other variables
    
    local(@bad_addresses) = ();
    $, = ', ';
    $" = ', ';
        
    # Open the connection to sendmail and line buffer output
    
    open(MAIL, "|$SENDMAIL -t");
    select(MAIL); 
    $| = 1; 
    select(STDOUT);
   
    # Output the message header
    
    print MAIL "To: @to\n";
    print MAIL "From: $from\n";
    print MAIL "CC: @cc\n" if $cc;
    print MAIL "BCC: @bcc\n" if $bcc;
    print MAIL "Subject: $subject\n";

    print MAIL "\n";

    # Output the message body.
        
    if ($body) {
        if (!($body =~ /^[\\\/:]/) && ($body =~ /\s/)) { print MAIL $body }
        elsif (-e $body && -T $body) { &parse_template($body, *MAIL) }
    }
    print MAIL "\n";

    # Close the connection
    
    close(MAIL);    
    return(0);
}


#--------------------------------------------------------------------------
# lock
#--------------------------------------------------------------------------
#
# Function:      Creates an exclusive lock for a file. The lock will       #
#                only work if other programs accessing the file are also   #
#                using this subroutine.                                    #
#                                                                          #
# Usage:         &lock($filename, $LOCK_DIR[, $MAX_WAIT]);                 #
#                                                                          #
# Variables:     $filename --   Name of file being locked.                 #
#                               Example "filename.html"                    #
#                $LOCK_DIR --   Path of directory to store lock files      #
#                               Should be "/tmp/" on UNIX sytems           #
#                               Example "/home/lockdir/"                   #
#                $MAX_WAIT --   Maximum seconds to wait if the file is     #
#                               already locked                             #
#                                                                          #
# Returns:       0 if successful                                           #
#                1 if $LOCK_DIR/$filename.tmp could not be created         #
#                2 if $filename is currently in use                        #
#                3 if lock file could not be created or opened             #
#                                                                          #
# Uses Globals:  $Error_Message for descriptive error messages             #
#                $NAME_LEN for maximum filename length                     #
#                                                                          #
# Files Created: Creates $LOCK_DIR/$filename.tmp                           #
#                Creates $LOCK_DIR/$filename.lok (exists only while file   #
#                  is locked)                                              #
#

sub lock {
    
    # Initialize variables
    
    local($filename, $LOCK_DIR, $MAX_WAIT) = @_; 
    local($wait, $lock_pid);
    local($temp_file) = "$LOCK_DIR$$.tmp";
    $Error_Message = '';
    
    local($lock_file) = $filename;
    $lock_file =~ tr/\/\\:.//d;         # Remove file separators/periods
    if ($NAME_LEN && ($NAME_LEN < length($lock_file))) {
        $lock_file = substr($lock_file, -$NAME_LEN);
    }
    $lock_file = "$LOCK_DIR$lock_file.lok";
    
    # Create temp file with PID
    
    if (!open(TEMP, ">$temp_file")) {
        $Error_Message = "Could not create $temp_file ($!).";
        return(1);
    }           
    print TEMP $$;
    close(TEMP);
    
    # Test for lock file
    
    if (-e $lock_file) {

        # Wait for unlock if lock file exists
        
        for ($wait = $MAX_WAIT; $wait; --$wait) {
            sleep(1);
            last unless -e $lock_file;
        }
    }
    
    # Check to see if there's still a valid lock
    
    if ((-e $lock_file) && (-M $lock_file < 0)) {
        
        # The file is still locked but has been modified since we started
                    
        unlink($temp_file);
        $Error_Message = "The file \"$filename\" is currently in use. Please try
again later.";
        return(2);
    }
    else {

        # There is either no lock or the lock has expired
        
        if (!rename($temp_file, $lock_file)) { 

            # Problem creating the lock file
        
            unlink($temp_file);
            $Error_Message = "Could not lock file \"$filename\" ($!).";
            return(3);
        }
        
        # Check to make sure the lock is ours

        if (!open(LOCK, "<$lock_file")) {
            $Error_Message = "Could not verify lock for file \"$filename\"
($!).";
            return(3);
        }
        $lock_pid = <LOCK>;
        close(LOCK);        
        if ($lock_pid ne $$) { 
            $Error_Message = "The file \"$filename\" is currently in use. Please
try again later.";
            return(2);
        }
        else { return(0) }
    }
}


#--------------------------------------------------------------------------
# unlock
#--------------------------------------------------------------------------
#
# Function:      Unlocks a file that has been locked using lock().         #
#                                                                          #
# Usage:         &unlock($filename, $LOCK_DIR);                            #
#                                                                          #
# Variables:     $filename --   Name of file being locked.                 #
#                               Example "filename.html"                    #
#                $LOCK_DIR --   Path of directory to store lock files      #
#                               Should be "/tmp/" on UNIX sytems           #
#                               Example "/home/lockdir/"                   #
#                                                                          #
# Returns:       0 if successful                                           #
#                1 if the lock file could not be deleted                   #
#                                                                          #
# Uses Globals:  $Error_Message for descriptive error messages             #
#                $NAME_LEN for maximum filename length                     #
#                                                                          #
# Files Created: Removes $LOCK_DIR/$filename.lok                           #
#

sub unlock {
    
    # Initialize variables
    
    local($filename, $LOCK_DIR) = @_;
    local($lock_file) = $filename;
    $Error_Message = '';
    
    $lock_file =~ tr/\/\\:.//d;         # Remove file separators/periods
    if ($NAME_LEN < length($lock_file)) {
        $lock_file = substr($lock_file, -$NAME_LEN);
    }
    $lock_file = "$LOCK_DIR$lock_file.lok";
    
    # Check to make sure the lock is ours

    if (!open(LOCK, "<$lock_file")) {
        $Error_Message = "Could not access the lock file for \"$filename\"
($!).";
        return(1);
    }
    $lock_pid = <LOCK>;
    close(LOCK);        
    if ($lock_pid ne $$) { 
        $Error_Message = "The file \"$filename\" is locked by another process.";
       return(2);
    }

    # Release the lock by unlinking the lock file
    
    if (!unlink($lock_file)) {
        $Error_Message = "Could not unlock file \"$filename\" ($!).";
        return(3);
    }
    return(0);
}


#--------------------------------------------------------------------------
# parse_template
#--------------------------------------------------------------------------
#
# Function:      Searches a file for variables identified as <<VARIABLE>>  #
#                and substitutes the value with the corresponding key from # 
#                the global associative arrays %VAR, %CONFIG, %FORM, or    #
#                %ENV (the arrays are checked in this order). Lines in the #
#                file beginning with '0:' will only be printed if a        #
#                variable appears on that line and has a value in one of   #
#                the four arrays.                                          #
#                                                                          #
# Usage:         &parse_template($template_file, *FILEHANDLE);             #
#                                                                          #
# Variables:     $template_file --  Full path to file containing template  # 
#                                   Example "/directory/file"              #
#                *FILEHANDLE --     Reference to filehandle to write       #
#                                   parsed file to.                        #
#                                   Example *FILE, *STDOUT                 #
#                                                                          #
# Returns:       0 if $template_file could not be opened                   #
#                1 if successful                                           #
#                                                                          #
# Uses Globals:  %VAR --            Miscellaneous variables assoc array    #
#                %CONFIG --         Configuration variables assoc array    #
#                %FORM --           Form variables assoc array             #
#                %ENV --            Environment variables assoc array      #
#                $Error_Message --  Set to text message if error           #
#

sub parse_template {
    local($template_file, *OUT) = @_;
    local($line, $line_copy, $changes);

    # Open the template file and parse each line
    
    if (!open(TEMPLATE, $template_file)) { 
        $Error_Message = "Could not open $template_file ($!).";
        return(0);
    }
    while ($line = <TEMPLATE>) {

        # Initialize our variables
        
        $line_copy = '';
        $changes = 0;
        
        # Search for variables in the current line
        
        while ($line =~ /<<([^>]+)>>/) {
            
            # Build up the new line with the section of $line prior to the 
            # variable and the value for $var_name (check %VAR, %CONFIG,
            # %FORM, then %ENV for match)
        
            ++$changes;
            if ($VAR{$1}) { $line_copy .= $` . $VAR{$1} }
            elsif ($CONFIG{$1}) { $line_copy .= $` . $CONFIG{$1} }
            elsif ($FORM{$1}) { $line_copy .= $` . $FORM{$1} }
            elsif ($ENV{$1}) { $line_copy .= $` . $ENV{$1} }
            else {
                --$changes;
                $line_copy .= $`;
            }
                 
            # Change $line to the section of $line after the variable
                
            $line = $';
        }
        
        # Set $line according to whether or not any matches were found
            
        $line = $line_copy ? $line_copy . $line : $line;
        
        # Print line depending on presence of 0: and variables existing
         
        if (($line !~ s/^0://) || !$line_copy || $changes) {
            print OUT $line;
        }
    }
    close(TEMPLATE);
    return(1);
}


#--------------------------------------------------------------------------
# format_date
#--------------------------------------------------------------------------
#
# Function:      Formats the date and time according to a format string.   #
#                                                                          #
# Usage:         &format_date($date[, $format, $gmt]);                     #
#                                                                          #
# Variables:     $date --       Date to format in time() format            #
#                               Example: 839120118                         #
#                $format --     Format specification string. The following #
#                               variables can be used within the string:   #
#                $gmt --        Set to non-null if GMT should be used      #
#                                                                          #
#                               <h> hour without padding (i.e. 6, 12)      #
#                               <0h> hour with padding (i.e. 06, 12)       #
#                               <mh> military hour w/ padding (i.e. 06,23) #
#                               <n> minutes without padding (i.e. 6, 35)   #
#                               <0n> minutes with padding (i.e. 06, 35)    #
#                               <s> seconds without padding (i.e. 6, 35)   #
#                               <0s> seconds with padding (i.e. 06, 35)    #
#                               <ap> am/pm according to hour               #
#                               <ap.> a.m./p.m. according to hour          #
#                               <time> time w/o seconds (i.e. 6:12 p.m.)   #
#                               <timesec> time w/secs (i.e. 6:12:14 p.m.)  #  
#                               <mtime> military time w/o secs (i.e. 18:12)#  
#                               <mtimesec> military time w/secs            #
#                                         (i.e. 18:12:14)                  #
#                               <m> month without padding (i.e. 6, 12)     #  
#                               <0m> month with padding (i.e. 06,12)       #
#                               <mon> 3-character month name (i.e. Jan)    #
#                               <month> full month name (i.e. January)     #
#                               <d> day of month w/o padding (i.e. 6, 23)  #
#                               <0d> day of month w/padding (i.e. 06, 23)  #
#                               <df> day of month with suffix (i.e. 23rd)  #
#                               <wd> weekday (0 - 7)                       #
#                               <wday> 3-character weekday name (i.e. Fri) #
#                               <weekday> full weekday name (i.e. Friday)  #
#                               <yr> 2-digit year (i.e. 96)                #
#                               <year> 4-digit year (i.e. 1996)            #
#                               <yd> day of year w/o padding (i.e. 6, 300) #
#                               <0yd> day of year w/padding (i.e. 006, 300)#
#                               <ydf> day of year with suffix (i.e. 300th) #
#                               <dst> dst if daylight savings time or blank#
#                               <dstnot> 'not ' or blank                   #
#                                                                          #
#                               In addition, writing a variable name in    #
#                               uppercase will result in any letters in    #
#                               the value of the variable to appear in     #
#                               uppercase.                                 #
#                               i.e. <WEEKDAY> will result in MONDAY       #
#                                                                          #
#                               Any text in $format that does not match a  #
#                               variable name will be left alone.          #
#                                                                          #
#                               Note: If <m> or <0m> appear after a : they #
#                               be interpreted as minutes, not month.      #
#                                                                          #
# Returns:       $format with all variables replaced w/appropriate values  #
#                If $format is omitted then an associative array will be   #
#                   returned with the lowercase variables (no brackets) as # 
#                   the keys and their corresponding values as the values. #  
#                A null string or array if $date is null or zero.          #
#                No other error checking is performed                      #
#                                                                          #
# Uses Globals:  @DAYS   -- Array containing full names of weekdays        #
#                           (Sunday first)                                 #
#                @MONTHS -- Array containing full names of months (January #
#                           first)                                         #
#

sub format_date {
    
    # Get the arguments
    
    local($date, $format, $gmt) = @_;
    local(@suffix) = ('th', 'st', 'nd', 'rd');
    local(%date_vars, $result, $var, $upper, $last_digit, $suffix);
    local($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);
    
    # Create default arrays if necessary
    
    if (!@DAYS) {
                @DAYS = ('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday',

                         'Friday', 'Saturday');
    }
    if (!@MONTHS) {
                @MONTHS = ('January', 'February', 'March', 'April', 'May',
'June',  
                           'July', 'August', 'September', 'October', 'November',
                          'December');
    }

    # Convert the date into local time or GMT if a date was given
    
    if (!$date) { return ($format ? '' : %date_vars) }
    if ($gmt) {
        ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
gmtime($date);
    }
    else {
        ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
localtime($date);
    }
    
    # Calculate all our variables
    
    $date_vars{'mh'} = sprintf("%02d", $hour);
    $date_vars{'mtime'} = sprintf("%02d:%02d", $hour, $min);
    $date_vars{'mtimesec'} = sprintf("%02d:%02d:%02d", $hour, $min, $sec);
    $date_vars{'ap'} = ($hour > 12) ? 'pm' : 'am';
    $date_vars{'ap.'} = ($hour > 12) ? 'p.m.' : 'a.m.';
    if ($hour > 12) { $hour -= 12 }
    elsif ($hour == 0) { $hour = 12 }
    $date_vars{'h'} = $hour;
    $date_vars{'0h'} = sprintf("%02d", $hour);
    $date_vars{'n'} = $min;
    $date_vars{'0n'} = sprintf("%02d", $min);
    $date_vars{'s'} = $sec;
    $date_vars{'time'} = sprintf("%d:%02d %s", $hour, $min, $date_vars{'ap.'});
    $date_vars{'timesec'} = sprintf("%d:%02d:%02d %s", $hour, $min, $sec,
$date_vars{'ap.'});
    $date_vars{'0s'} = sprintf("%02d", $sec);
    $date_vars{'m'} = $mon + 1;
    $date_vars{'0m'} = sprintf("%02d", $mon + 1);
    $date_vars{'mon'} = substr($MONTHS[$mon], 0, 3);
    $date_vars{'month'} = $MONTHS[$mon];
    $date_vars{'d'} = $mday;
    $date_vars{'0d'} = sprintf("%02d", $mday);
    $last_digit = substr($mday, -1);
    if ((($mday < 11) || ($mday > 20)) && $last_digit < 4) {
        $date_vars{'df'} = $mday . $suffix[$last_digit];
    }
    else { $date_vars{'df'} = $mday . 'th' }
    $date_vars{'wd'} = $wday;
    $date_vars{'wday'} = substr($DAYS[$wday], 0, 3);
    $date_vars{'weekday'} = $DAYS[$wday];

# These calculations do not work for 2000.
#    $year = sprintf("%02d", $year);
#    $date_vars{'yr'} = $year % 100;
# Try this:
     $year = 1900 + $year;
     $date_vars{'yr'} = substr($year,2,2);

    $date_vars{'year'} = 1900 + $year;
    $date_vars{'yd'} = $yday;
    $date_vars{'0yd'} = sprintf("%03d", $yday);
    $last_digit = substr($yday, -1);
    if ((($yday % 100 < 11) || ($yday % 100 > 20)) && $last_digit < 4) {
        $date_vars{'ydf'} = $yday . $suffix[$last_digit];
    }
    else { $date_vars{'ydf'} = $yday . 'th' }
    $date_vars{'dst'} = $isdst ? 'dst' : '';
    $date_vars{'dstnot'} = $isdst ? '' : 'not ';

    # Scan the format string for variables, replacing them as we go
    
    while ($format =~ /<([^>]+)>/) {
    
        # Check if the variable name is in uppercase
        
        $_ = $1;
        $upper = tr/A-Z/a-z/;
        
        # Make an effort to correct m where n was intended

        if (substr($_, -1) eq 'm') {
            if ($` eq ':') { substr($_, -1) = 'n' };
        }

        # If the variable name is in uppercase convert the value to uppercase
        
        if ($upper) { 
            $var = $date_vars{$_};
            $var =~ tr/a-z/A-Z/;
            $result .= $` . $var;
        }
        else { $result .= $` . $date_vars{$_} }
        
        # Change format to the rest of the string after the match
        
        $format = $';
    }

    # Set $result according to whether or not any matches were found
    
    $result .= $format;
    
    # Return the result or the variable array depending how the routine was called
    
    if ($result) { $result }
    else { %date_vars }
}
