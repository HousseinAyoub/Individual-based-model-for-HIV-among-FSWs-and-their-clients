clust = parcluster();
job = createCommunicatingJob(clust,'Type','Pool');
task = createTask(job, @Multiple_Runs_for_Average, 79, {1}); %112 is the number of ourputs (HIAM please check this) and {1} is the number of estimated parameter  
tic; submit(job); toc
