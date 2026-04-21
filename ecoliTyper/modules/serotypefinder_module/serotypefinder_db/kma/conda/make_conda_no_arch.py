from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap
from io import StringIO

# Create an instance of YAML
yaml = YAML()
yaml.default_flow_style = False  # Set the default flow style to block style

# Create an ordered dictionary for each section
package = CommentedMap()
package['name'] = 'kma'

# Read the version from version.h
with open('version.h', 'r') as f:
    for line in f:
        if line.startswith('#define'):
            package['version'] = line.split()[2].replace("\"", "")

source = CommentedMap()
source['url'] = 'https://bitbucket.org/genomicepidemiology/{}/get/{}.tar.gz'.format(package['name'], package['version'])

build = CommentedMap()
#build['number'] = 0
build['noarch'] = 'generic'

requirements = CommentedMap()
requirements['build'] = ['make', '{{ compiler("c") }}']
requirements['host'] = ['zlib']
requirements['run'] = ['zlib']

about = CommentedMap()
about['home'] = 'https://bitbucket.org/genomicepidemiology/kma'
about['summary'] = 'KMA is mapping a method designed to map raw reads directly against redundant databases, in an ultra-fast manner using seed and extend.'
about['license'] = 'Apache-2.0'

extra = CommentedMap()
identifiers = CommentedMap()
identifiers['doi'] = '10.1186/s12859-018-2336-6'
extra['identifiers'] = identifiers

# Create a dictionary for the entire YAML content
data = CommentedMap()
data['package'] = package
data['source'] = source
data['build'] = build
data['requirements'] = requirements
data['about'] = about
data['extra'] = extra

# Use StringIO to capture the YAML string
stream = StringIO()
yaml.dump(data, stream)
yaml_str = stream.getvalue().replace("\"{{", "{{").replace("}}\"", "}}")

# Print the YAML string
print(yaml_str)

# Write the YAML string to a file
with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
