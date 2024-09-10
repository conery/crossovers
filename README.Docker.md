# Docker Image

The crossover explorer scripts in this repo are available in a Docker image.

We do not recommend using Docker to run the `xo peaks` command to create the set of blocks.  In our initial tests the command runs very slowly and consumes an huge amount of the host system's resources.

The `xo view` command, however, runs well, and is very easy to set up and run. 

## Overview

If you're not familiar with Docker:

- A Docker **image** is basically a virtual machine, a "snapshot" that contains a scaled-down Linux operating system (with only enough functionality to run Python programs), a Python interpreter, the library modules needed by our program, and our Python code
- A Docker **container** is created when we want to run the code in an image.  A container is based on an image, and includes all the computing resources needed to run the code.  In our case, the resources we need are a directory that has the data files and a network port that allows us to view the GUI produced by the program running in the container.

In the steps below you'll download ("pull") the image that contains our application, and then you'll create a container based on the image.  The step that creates the container also starts it.

When you're done using the container you can stop it.  The next time you want to view the data you should **restart** the container -- if you repeat the steps that create a container you'll end up with two containers and it's easy to get confused.

#### A Note About Docker Accounts

You do not need to log in to Docker in order to pull an image.  Simply start Docker and use the process described below to download the image for this project.  If you have a Docker account you can log in, or you can log in using the lab account, but that is only necessary if you want to manage the images that are stored on Docker Hub.

## Pull the Image

Open Docker Desktop by clicking on the whale icon in your toolbar and select "Go to the Dashboard".

> If the whale icon is not there it means Docker is not running.  Start Docker by double-clicking on the Docker application in your Applications folder.

Find the search box in the menu bar at the top of the window.  Type "libudalab" and hit return.  You should see a list of names of Libuda Lab images.

Select `libudalab/crossovers` from the table.  The window will show a box labeled "Tag".  Make sure "latest" is selected, then click Pull.

#### A Note About Image Versions

In the future, if there is an update to the image, you can repeat the steps above to download the newest version.  If for some reason you need an older version you can select it from the menu instead of "latest".

## Data

The container you're going to create in the next step expects to find the data files it uses in to be in a single folder.   The default names for the files are:

- `peaks.csv`:  a CSV file with descriptions of blocks of SNPs found by the `xo peaks` script

- `BSP_TIGER.intervals_dataframe.pickle.gzip`:  a summary of the regions of each chromosome, a Pandas data frame where lines describe regions of consecutive SNPs

If you want to use different file names you need to specify those names when you create the container, described in the next section.  You can use any names you like, as long as (a) the peaks data is in a plain text CSV file, (b) the interval data is a "pickled" data frame compressed with GZIP, and (c) both files are in the same folder.

## Create the Container

In the Docker Dashboard window click on the Images tab in the toolbar on the left side of the window, and make sure Local is selected (instead of Hub).

You'll see a list of images you have pulled from DockerHub.  Find the line for `libudalab/crossovers` and click the run button (right-facing triangle) on the right side of the line.

A window will pop up to allow you to specify options.  Click the downward pointing triangle to show the options settings.

- Click in the box that says "container name" a give your new container a name, something like "crossovers"

- Under Ports, click in the space labeled "host port" and enter 8000
- Under Volumes, click in the space labeled "host path".  A file browser will pop up.  Navigate to the folder where you stored the data files and click Open
- Also under Volumes, click the space labeled "container path" and enter `/data`
- (Optional)  If you have different names for your data files than the ones shown above you need to enter those names in the section called  Environment variables.  In the space labeled Variable enter `XO_PEAKS`, and in the space labeled Value enter the name of the CSV file that has your peaks data.  Then click the plus button to add another set of Variable and Value boxes.  Type `XO_INTERVALS` in the variable box and the name of your intervals file in the value box.

Click the Run button to start the container.

### Container View

When Docker starts the container the display in the Docker Desktop changes to show you the output printed on the terminal by the application running in the container.  You'll see it printing status messages (like "reading data") and eventually a message that says "Launching server at http://localhost:5006".

## Connect to the GUI

Open a browser window and enter this URL:  `http://localhost:8000` .

## Closing the GUI

When you're done:

- Stop the container:  Go to the Docker Dashboard and click on Containers in the sidebar.  You should see your running container with the name you gave it.  On the right side of the line there is a square stop button -- click that button to stop the container.
- Close the browser window.

## Restarting the Container

The next time you want to look at the data you should restart the container.  Simply go to the Docker Desktop and click on Containers on the left side of the window.  Find your container and look for the start button on the right side of the line (after you start the container the button will change into a stop button; come back here and click that button when you're done).

## Troubleshooting

Here are some things that can go wrong.

### Port Conflict

If you're using port 8000 for some other application, _e.g._ if you have a Jupyter notebook open, Docker will display an error when it tries to create the container.  Simply choose a different port number, like 8001, and enter it in the options form.  Remember to use this port number in your browser, _e.g._ you'll go to `localhost:8001`.

### To Be Continued...

