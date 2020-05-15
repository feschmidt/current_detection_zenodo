# %%
import pkg_resources

installed_packages = [(d.project_name, d.version)
                      for d in pkg_resources.working_set]
# %%
installed_packages = sorted(installed_packages, key=lambda tup: tup[0].lower())

# %%
with open('python_info.txt', 'w') as f:
    for item in installed_packages:
        f.write("{}=={}\n".format(item[0], item[1]))
