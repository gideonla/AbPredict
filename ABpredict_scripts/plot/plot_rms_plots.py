import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import pdb
from matplotlib.legend_handler import HandlerLine2D



filename=sys.argv[1]
top_cluster=sys.argv[2]

with open(top_cluster) as f:
    content = f.readlines()
content = [x.strip() for x in content] 
print "content",content

#f = open(filename, 'r')

data = pd.read_csv(filename, sep=" ")#, index="dates", names=["energy", "whole_vl_vh", "L1", "L2", "L3", "H1", "H2","H3"])
print data
top_cluster_list=[]
i=0
while i < len(data.ix[:,0]):
	print(i,data.ix[i,0])
	for j in content:
		#print j
		#print (data.ix[i,0])
		if j in  data.ix[i,0]:
			#print j
			#print data.ix[i,0]
			top_cluster_list.append(i)
	i=i+1
print (top_cluster_list)
print (min(data.ix[:,1]))

# make plot figure of energy vs. whole model

colors = (0,0,0)


def make_plot(num, name):
	f = plt.figure(num, figsize=(10,10))
	plt.scatter(data.ix[:,num], data.ix[:,1],s=100, c=colors, alpha=0.5,label='_nolegend_') #plot x-axis = whole rmsd vs. y-axis = energy
	model1=plt.scatter(data.ix[top_cluster_list[0],num], data.ix[top_cluster_list[0],1],color="red", alpha=1,s=150, label="model1")
	model2=plt.scatter(data.ix[top_cluster_list[1],num], data.ix[top_cluster_list[1],1],color="orange", alpha=1,s=150, label="model2")
	model3=plt.scatter(data.ix[top_cluster_list[2],num], data.ix[top_cluster_list[2],1],color="green", alpha=1,s=150, label="model3")
	plt.title('Energy vs. %s CDR C=O rmsd'%(name.split("_")[0]),fontsize=20)
	plt.xlim([min(data.ix[:,num])-0.1,max(data.ix[:,3])+0.1])
	plt.ylim([min(data.ix[:,1])-1,max(data.ix[:,1])+10])
	plt.xlabel('RMSD', labelpad=-1,fontsize=20)
	plt.ylabel('R.e.u',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=20)
	plt.tight_layout()
	plt.legend(loc='upper left',markerscale=0.5, scatterpoints=1)

	plt.savefig(name+"_vs_reu.png")		

graph_titles = ["whole_rmsd","L1_rmsd","L2_rmsd","L3_rmsd","H1_rmsd","H2_rmsd","H3_rmsd"]
for i in range(2,len(graph_titles)+2):
	make_plot(i,graph_titles[i-2])

#plot tile angle
f = plt.figure(10, figsize=(10,10))
angle_mean=data.ix[:,9].mean()
new_angle = data.ix[:,9] - angle_mean
plt.scatter(new_angle[:], data.ix[:,1],s=100, c=colors, alpha=0.5,label='_nolegend_') #plot x-axis = whole rmsd vs. y-axis = energy
plt.title('Energy vs. vl/vh angle',fontsize=20)
#plt.set_xlim([min(data.ix[:,9])-0.5,max(data.ix[:,9])+0.5])
#plt.set_ylim([min(data.ix[:,1])-50,max(data.ix[:,1])+50])
plt.xlabel('Angle', labelpad=-1,fontsize=20)
plt.ylabel('R.e.u',fontsize=20)

plt.scatter(new_angle[top_cluster_list[0]], data.ix[top_cluster_list[0],1],color="red", alpha=1,s=150, label="model1")
plt.scatter(new_angle[top_cluster_list[1]], data.ix[top_cluster_list[1],1],color="orange", alpha=1,s=150, label="model2")
plt.scatter(new_angle[top_cluster_list[2]], data.ix[top_cluster_list[2],1],color="green", alpha=1,s=150, label="model3")
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()
plt.legend(loc='upper left',markerscale=0.5, scatterpoints=1)
plt.savefig("tilt_vs_reu.png")		





