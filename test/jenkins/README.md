# Jenkins Setup Instructions and Usage Information

This file describes SUNDIALS [Jenkins CI](https://www.jenkins.io/) setup.

## Installing Jenkins

The following steps are adapted from the Jenkins Red Hat installation
[instructions](https://www.jenkins.io/doc/book/installing/linux/).

1. Install Java 11

   ```
   sudo yum install java-11-openjdk-devel
   ```

2. Install the stable version

   ```
   sudo wget -O /etc/yum.repos.d/jenkins.repo https://pkg.jenkins.io/redhat-stable/jenkins.repo
   sudo rpm --import https://pkg.jenkins.io/redhat-stable/jenkins.io.key
   sudo yum upgrade
   sudo yum install jenkins
   sudo systemctl daemon-reload
   ```

   or weekly version of Jenkins

   ```
   sudo wget -O /etc/yum.repos.d/jenkins.repo https://pkg.jenkins.io/redhat/jenkins.repo
   sudo rpm --import https://pkg.jenkins.io/redhat/jenkins.io.key
   sudo yum upgrade
   sudo yum install jenkins
   sudo systemctl daemon-reload
   ```

3. Start Jenkins

   ```
   sudo systemctl start jenkins
   ```

4. Check Jenkins status

   ```
   sudo systemctl status jenkins
   ```

5. Open a web browser and navigate to `http://system_address:8080/` where
   `system_address` is the network address of the machine Jenkins was installed
   on or `localhost`.

6. Click `Install suggested plugins`.

7. When the installation is complete, fill in the information needed
   to create the first admin user then click `Save and Continue`.

8. Edit the jenkins URL if desired and click `Save and Continue`.

9. Click `Start using Jenkins`.

## Additional Plugins

To install plugins:

1. Navigate to Manage Jenkins > Manage Plugins

2. Click the "Available" tab

3. Search for the plugin name

4. Check the box next to plugin

5. Click "install without restart"

Suggested plugins:

* Bitbucket Branch Source -- provides better integration with Bitbucket server,
  Git branches, and pull requests.

* Collapsing Console Sections -- provides better readability of the
  Jenkins console output by enabling collapsible sections.

* Dashboard view -- Customizable dashboard that can present various views of job
  information.

* Rebuilder -- allows the user to rebuild a parameterized build without entering
  the parameters again.



## Setup Collapsing Console Sections

The Collapsing Console Sections Plugin works by parsing the Jenkins console
output to find lines, defined by Java regular expressions, that mark the start
and end of a collapsible section.

1. Navigate to Manage Jenkins > Configure System

2. Under "Collapsing Console Sections" check the box next to Enable
   sections numbering.

3. To create a Jenkins Setup section, click on "Add Console Section"
   and enter/do the following:

   Section name:        Jenkins Setup
   Section starts with: ((Branch indexing)|(Started by user .*))
   Section ends with:   END JENKINS SETUP
   Check the box next to "Collapse Sections by default"

4. To create a CMAKE section, click on "Add Console Section" and
   enter/do the following:

   Section name:        CMAKE
   Section starts with: START CMAKE
   Section ends with:   (cmake returned .*)
   Check the box next to "Collapse Sections by default"

5. To create a BUILDING section, click on "Add Console Section" and
   enter/do the following:

   Section name:        BUILDING
   Section starts with: START MAKE
   Section ends with:   (make returned .*)
   Check the box next to "Collapse Sections by default"

6. To create a TESTING section, click on "Add Console Section" and
   enter/do the following:

   Section name:        {1}
   Section starts with: TEST: (.*)
   Section ends with:   (PASSED|FAILED): (.*)
   Check the box next to "Collapse Sections by default"

7. To create a INSTALLING section, click on "Add Console Section" and
   enter/do the following:

   Section name:        INSTALLING
   Section starts with: START INSTALL
   Section ends with:   (make install returned .*)
   Check the box next to "Collapse Sections by default"

9. Click "Save"

## Adding Credentials

For Jenkins to use the Bitbucket server API we need to provide it with login
credentials.

1. Navigate to Manage Jenkins > Manage Credentials > Jenkins > Global credential

2. Click "Add Credentials"

3. The Kind should be "Username and password" and the Scope should be
   "Global (Jenkins, nodes, items, all child items, etc)". Fill in the username
   and password then click OK.

## Adding a Bitbucket Server:

Because the repository is located in a local Bitbucket server we need to
give this information to Jenkins:

1. Navigate to Manage Jenkins > Configure System

2. Under Bitbucket Endpoints, click on "Add" and select "Bitbucket Server"

3. Enter a name for the server and the server URL

4. Click "Save"

## Selecting the Git Installation

To change which installation of Git Jenkins uses

1. Navigate to Manage Jenkins > Global Tool Configuration

2. Under Git, change "Path to Git executable" to the desired installation

3. Click "Save"

NOTE: This setting only effects what version of Git Jenkins will used to
interact with the repo. If any regression test scripts use Git commands they
will use the default install of Git unless a newer version is added to the
path or called explicitly. See below for updating the `PATH` environment
variable.

## Updating Environment Variables

If regression testing scripts require that certain path are set in the `PATH`
environment variable, the following steps detail how to modify
PATH for each Jenkins testing node.

1. Navigate to Manage Jenkins > Manage Nodes

2. Click the gear icon next to the node you want to configure

3. Under "Node Properties" check the box next to "Environment variables"

4. Click "Add" and enter the PATH+SOMENAME in the Name field and the desired
   path in the Value field. For example to add the path for a newer version of
   git to PATH do enter the following:

   Name:  PATH+NFS_GIT
   Value: /usr/apps/git/2.9.4/bin/git

5. Click "Save"

## Setting the Number of Build Executors

Depending on the machine Jenkins is installed on it may or may not be desirable
to allow multiple regression tests to run in parallel on a node (machine). The
number of jobs Jenkins can run in parallel can be adjusted on a per node basis.
If there are more jobs to run than idle nodes Jenkins will place the jobs in a
first in first out queue and run them as resources become available.

1. Navigate to Manage Jenkins > Manage Nodes

2. Click the gear icon next to the desired node

3. Change "# of executors" to the desired number of test jobs that can run in
   parallel

4. Click "Save"

## Enabling Anonymous Read Access

For users to view the build status without needing to logging in to Jenkins will
enable anonymous read access.

1) Navigate to Manage Jenkins > Configure Global Security

2) Under Access Control > Authorization, check the box in front of
   "Allow anonymous read access"

3) Click "Save"

## Changing the Jenkins User

By default Jenkins is installed and runs as the local user jenkins. In some
cases is the desirabel to have Jenkins run as a different user. The following
instructions describe how to change the Jenkins user.

1. Stop Jenkins

   ```
   sudo systemctl stop jenkins
   ```

2. Open the Jenkins configuration file, change the `JENKINS_USER` variable to
   the desired user name, and save the revised file

   ```
   sudo nano /etc/sysconfig/jenkins
   JENKINS_USER="username"
   ```

3. Change the ownership of the Jenkins home, webroot, and logs

   ```
   sudo chown -R username:username /var/lib/jenkins
   sudo chown -R username:username /var/cache/jenkins
   sudo chown -R username:username /var/log/jenkins
   ```

4. Start Jenkins

   ```
   sudo systemctl start jenkins
   ```

## Adding a Jenkins Pipeline

A Jenkins Pipeline is domain specific language for implementing a project's
entire build/test/deploy pipeline in a Jenkinsfile that is stored in the
project's code repo. The following steps setup Jenkins pipeline for a
multibranch Bitbucket server repository.

1. Click on "New Item"

2. Enter an item name e.g. SUNDIALS

3. Click "Multibranch Pipeline" then click "OK"

4. Under Branch Source, click "Add source" and select "Bitbucket"

5. Under Bitbucket set the following:

   a. Server: select the server name entered when adding the Bitbucket
      server above

   b. Credentials: select the user name entered for accessing the
      Bitbucket server API

   c. Owner: sundials

   d. Repository Name: select the desired repo name (e.g. sunrepo)

6. Under "Build Configuration" provide the path to the Jenkinsfile

7. Click "Save"

Note: If the hooks to notify Jenkins of new pushes to the Bitbucket are not
working, under" Scan Multibranch Pipeline Trigger" check the box next to
"Periodically if not otherwise run" and set the desired interval.

Jenkins will now scan the repository for branches and pull requests that contain
a Jenkinsfile and run the pipeline for any branches or pull requests with a
Jenkinsfile.
