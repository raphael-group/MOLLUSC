import treeswift



#input: file name
file_name = "testing.txt"
f = open(file_name, 'r')
data = f.read()
newick_tree_start = data.find('[&R] ') + 5
new_tree_end = data.find(';', newick_tree_start)
tree_string = data[newick_tree_start:new_tree_end+1]
# print(tree_string)

newick_tree = treeswift.read_tree(tree_string, 'newick')
print(newick_tree.height())
newick_tree.draw(show_labels=True)



