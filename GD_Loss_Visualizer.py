from __init__ import *
import matplotlib.pyplot as plt
import re 

def extract_numbers_in_parentheses(string):
    pattern = r"\((\d+)\)"
    matches = re.findall(pattern, string)
    return matches

def get_gd_node_to_species_loss_dic(sp_dic):
    loss_dic = {'No-duplicate loss': [], 'One-duplicate loss': [], 'Two-duplicate loss': []}
    for k,v in sp_dic.items() :
        s=extract_numbers_in_parentheses(k)
        s_set=set(s)
        if s_set=={'2'}:
            loss_dic['No-duplicate loss'].append(v)
        elif s_set=={'2','1'} or s_set=={'2','1','0'} or s_set=={'1','0'}:
            loss_dic['One-duplicate loss'].append(v)
        else:
            loss_dic['Two-duplicate loss'].append(v)

    sum_loss_dic = {key: sum(value) for key, value in loss_dic.items()}

    return sum_loss_dic


def visualizer_sum_loss_dic(sum_loss_dic, sps, gd_id, out_dir):
    keys = list(sum_loss_dic.keys())
    values = list(sum_loss_dic.values())
    color = 'lightblue'

    plt.figure(figsize=(10, 6))  

    bars = plt.bar(keys, values, color=color)

    plt.ylabel('GD num', fontsize=14) 
    plt.title('Count of ' + sps + ' Species under Node ' + gd_id, fontsize=16)  
    plt.yticks(fontsize=12)  

    # 添加数据标签
    for key, value in zip(keys, values):
        plt.text(key, value, str(value), ha='center', va='bottom', fontsize=12)

    plt.tight_layout()  
    plt.savefig(f'{out_dir}/{gd_id}_{sps}.png')
    plt.cla()
    plt.close("all")  

def generate_plt(input_folder,out_dir):

	input_dirs=os.listdir(input_folder)
	pbar = tqdm(total=len(input_dirs), desc="Processing file", unit="file")
	for file in input_dirs:
		pbar.set_description(f"Processing {file}")
		full_path = os.path.join(input_folder, file)
		new_dic=read_and_return_dict(full_path)
		sum_loss_dic=get_gd_node_to_species_loss_dic(new_dic)
		gd_id,sps=file.split('.')[0].split('_')
		visualizer_sum_loss_dic(sum_loss_dic,sps,gd_id,out_dir)
		pbar.update(1)
	pbar.close()
if __name__ == "__main__":
    out='outfile'
    input_folder='input'
    os.makedirs(out, exist_ok=True)
    enerate_plt(input_folder,out)

