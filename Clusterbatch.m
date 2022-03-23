clust = parcluster();
job = createCommunicatingJob(clust,'Type','Pool');
task = createTask(job, @Multiple_Runs_for_Average, 108, {2}); %112 is the number of ourputs and {2} is the number of estimated parameters 
tic; submit(job); toc
 