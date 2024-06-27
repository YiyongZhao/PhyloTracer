from __init__ import *

def process_gd_result(gd_file):
	gds=[]
	with open(gd_file,'r') as f:
		for i in f :
			if i.startswith('#'):
				continue
			a=i.strip().split('\t')
			gds.append((a[1],a[5]))
	return gds

def get_count_dic(gds):
	gd_set=set()
	count_dic={}
	for i in gds:
		gd_id,level=i
		if gd_id in gd_set:
			continue
		else:
			if level in count_dic:
				count_dic[level]+=1
			else:
				count_dic[level]=1
			gd_set.add(gd_id)
	return count_dic

def mark_sptree(sptree:object,count_dic:dict)->object:
	num_tre_node(sptree)
	sptree.ladderize()
	sptree.sort_descendants("support")
	for node in sptree.traverse():

		nstyle = NodeStyle()
		nstyle["fgcolor"] = "black"
		nstyle["size"] = 0
		nstyle["shape"] = "circle"
		node.set_style(nstyle)

		num=str(count_dic.get(node.name,0))
		node.add_face(TextFace(num+' GD', fsize=5, fgcolor="red"), column=0, position="branch-top")
	return sptree.render('phylotracer_gd_visualizer.PDF')

def gd_visualizer_main(sptree,gd_result):
	gds=process_gd_result(gd_result)
	count_dic=get_count_dic(gds)
	mark_sptree(sptree,count_dic)



