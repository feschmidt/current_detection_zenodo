# outputs directory structure excluding certain file extensions to tree_info
# run via `bash tree_info.sh`
tree -I '*.pkl|*.png|*.pyc|__pycache__|__init__.py' > tree_info.txt
