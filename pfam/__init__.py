import re
import requests
import xml.etree.ElementTree as ElementTree
from datetime import datetime


def request (path: str, params: dict = {}, headers: dict = None, timeout: float = None) -> ElementTree:

    # Default value for option 'output'
    params.setdefault('output', 'xml')

    # Compose the url
    url = 'https://pfam.xfam.org' + path 

    # Make the request to pfam
    response = requests.get(url, params=params, headers=headers, timeout=timeout)

    # Throw an execption if HTTP status is not ok
    response.raise_for_status ()

    # If the output is an XML file, parse it
    if params['output'] == 'xml':

        # Remove the namespace from the xml
        xml = re.sub(' xmlns="[^"]+"', '', response.text, count=1)

        # Parse the XML and create the root element
        root = ElementTree.fromstring(xml)

        # Check if is an error
        if root.tag == 'error':
            raise Exception(root.text)

        # Process all values in the XML tree
        xml_process_tree(root)

        # Return the XML root
        return root

    # Return the response as text
    else:
        return response.text


def xml_process_value(value: str):

    # The value can be None
    if value is None:
        return None

    # This function process only strings 
    if type(value) != str:
        return value

    # Otherwise it is a string
    value = value.strip()

    # Trasform into None if is the empty string
    if value == '':
        return None

    # Try to convert into a float
    try:
        return float(value)

    except ValueError:
        pass

    # Return the string stripped
    return value


def xml_process_tree(element: ElementTree):

    # Visit all element in the XML tree
    for el in element.iter():

        # Process the text contained in the element
        el.text = xml_process_value(el.text)

        # Process all the attributes
        for name in el.attrib:
            el.attrib[name] = xml_process_value(el.attrib[name])

class PfamEntry:

    def __init__(self, type: str, id: str, accession: str, description: str):

        # Set the members as the parameters
        self.id          : str = id
        self.accession   : str = accession
        self.description : str = description
        self.type        : str = type

    def fetch (self):

        if self.type == 'Pfam-A':
            return family(self.accession)

        elif self.type == 'Clan':
            return clan(self.accession)

        elif self.type == 'sequence':
            return protein(self.accession)

class PfarmGOTerm:

    def __init__(self, id: str, category: str, text: str):

        # Save the parameters as members
        self.id       : str = id
        self.category : str = category
        self.text     : str = text
        

class PfamCurationDetails:

    def __init__(self, element: ElementTree):

        # Typed members of this object
        self.status              : str = None
        self.seed_source         : str = None
        self.num_archs           : int = None
        self.num_species         : int = None
        self.num_structures      : int = None
        self.percentage_identity : float = None
        self.av_length           : float = None
        self.av_coverage         : float = None
        self.type                : str = None

        # Process child nodes
        for el in element:
        
            # A special parsing for <num_seqs> tag
            if el.tag == 'num_seqs':
                self.num_seqs_seed: float = el.find('seed').text
                self.num_seqs_full: float = el.find('full').text

            # Set an attribute of this object using the text content of the element as value
            else:
                setattr(self, el.tag, el.text)


class PfamCutoff:

    def __init__(self, sequence: float, domain: float):

        # Save the parameters as members
        self.sequence : float = sequence
        self.domain   : float = domain


class PfamHmmDetails:

    def __init__(self, element: ElementTree):

        # Attributes of the <hmm_details> tag
        self.hmmer_version : str = element.get('hmmer_version')
        self.model_version : str = element.get('model_version')
        self.model_length  : int = element.get('model_length')

        # Other propperties fetched form the child nodes
        self.build_commands  : str  = None
        self.search_commands : str  = None
        self.cutoffs         : dict = { }

        # Process child nodes
        for el in element:

            # A special parsing for <cutoffs> tag
            if el.tag == 'cutoffs':
                for cutoff in el:
                    self.cutoffs[cutoff.tag] = PfamCutoff(cutoff.find('sequence').text, cutoff.find('domain').text)

            # Set an attribute of this object using the text content of the element as value
            else:
                setattr(self, el.tag, el.text)
        

class PfamFamily:

    def __init__(self, release_version: str, release_date: datetime, entry: ElementTree):

        # The release version of the pfam database
        self.pfam_release_version : str      = release_version
        self.pfam_release_date    : datetime = release_date

        # Text of the family
        self.description : str = None
        self.comment     : str = None
        
        # Process child nodes
        for el in entry:
            
            # Trasform the element curation_details in a object
            if el.tag == 'curation_details':
                self.curation_details: PfamCurationDetails = PfamCurationDetails(el)

            # Trasform the element hmm_details in a object
            elif el.tag == 'hmm_details':
                self.hmm_details: PfamHmmDetails = PfamHmmDetails(el)
            
            elif el.tag == 'clan_membership':
                self.clan : PfamEntry = PfamEntry('Clan', el.get('clan_id'), el.get('clan_acc'), None)

            # Parse GO terms in a list of dicts
            elif el.tag == 'go_terms':
                self.go_terms: [PfarmGOTerm] = []

                # Iterate over categories
                for cat in el:

                    # The category name
                    category = cat.get('name')

                    # Insert in the list all the terms
                    for term in cat:
                        self.go_terms.append(PfarmGOTerm(term.get('go_id'), category, term.text))
            
            # Set an attribute of this object using the text content of the element as value
            else:
                setattr(self, el.tag, el.text)

        # Entry of this family (for consistency)
        self.entry: PfamEntry = PfamEntry(entry.get('entry_type'), entry.get('id'), entry.get('accession'), self.description)

class PfamMatch:

    def __init__(self, element: ElementTree):

        # Attribute data of the tag <matches>
        self.entry: PfamEntry = PfamEntry(element.get('type'), element.get('id'), element.get('accession'), None)

        # List of locations of matches
        self.locations: [PfamLocation] = []

        # Populate the list of locations
        for loc in element:
            self.locations.append (PfamLocation(loc))

    def fetch_family (self) -> PfamFamily:

        # Request the family to the server
        return self.entry.fetch()
            


class PfamLocation:

    def __init__(self, element: ElementTree):

        # Set members using attributes of <location> tag
        self.start       : int = element.get('start')
        self.end         : int = element.get('end')
        self.ali_start   : int = element.get('ali_start')
        self.ali_end     : int = element.get('ali_end')
        self.hmm_start   : int = element.get('hmm_start')
        self.hmm_end     : int = element.get('hmm_end')
        self.bitscore    : float = element.get('bitscore')
        self.evalue      : float = element.get('evalue')
        self.evidence    : str = element.get('evidence')
        self.significant : str = element.get('significant')

        self.hmm          : str = None
        self.match_string : str = None
        self.pp           : str = None
        self.seq          : str = None
        self.raw          : str = None

        # Process child nodes
        for el in element:

            # Set an attribute of this object using the text content of the element as value
            setattr(self, el.tag, el.text)

class PfamProtein:

    def __init__(self, release_version: str, release_date: datetime, entry: ElementTree):

        # The release version of the pfam database
        self.pfam_release_version : str      = release_version
        self.pfam_release_date    : datetime = release_date

        # Attributes of this entry
        self.db_name            : str   = entry.get('db')
        self.db_release_version : float = entry.get('db_release')

        # Text of the family
        self.description : str = None
        self.comment     : str = None

        self.taxonomy_id  : str   = None
        self.species_name : str   = None
        self.taxonomy     : [str] = None
        self.sequence     : str   = None

        # List of matches
        self.matches: [PfamMatch] = []

        # Process child nodes
        for el in entry:

            if el.tag == 'taxonomy':
                self.taxonomy_id  : int   = el.get('tax_id')
                self.species_name : str   = el.get('species_name')
                self.taxonomy     : [str] = el.text.rstrip('.').split('; ')

            elif el.tag == 'sequence':
                self.sequence_version : int = el.get('version')
                self.sequence         : str = el.text
                
            elif el.tag == 'matches':
                for match in el:
                    self.matches.append(PfamMatch(match))
            
             # Set an attribute of this object using the text content of the element as value
            else:
                setattr(self, el.tag, el.text)

        # Entry of this protein
        self.entry: PfamEntry = PfamEntry(entry.get('entry_type'), entry.get('id'), entry.get('accession'), self.description)

    def fetch_family (self) -> PfamFamily:

        # Only if it has a match is possible
        if len(self.matches) > 0:
            return self.matches[0].fetch_family()

class PfamClanMember:
    
    def __init__(self, element: ElementTree):

        self.entry           : PfamEntry = PfamEntry('Pfam-A', element.get('id'), element.get('accession'), None)
        self.num_occurrences : float     = element.get('num_occurrences')
        self.percentage_hits : float     = element.get('percentage_hits')

class PfamClan:

    def __init__(self, release_version: str, release_date: datetime, entry: ElementTree):

        # The release version of the pfam database
        self.pfam_release_version : str      = release_version
        self.pfam_release_date    : datetime = release_date

        # Text of the family
        self.description : str = None
        self.comment     : str = None

        # Members
        self.families = []
        
        # Process child nodes
        for el in entry:
            
            # Trasform the element curation_details in a object
            if el.tag == 'members':
                for member in el:
                    self.families.append (PfamClanMember(member))
            
            # Set an attribute of this object using the text content of the element as value
            else:
                setattr(self, el.tag, el.text)

        # Entry of this clan
        self.entry: PfamEntry = PfamEntry(entry.get('entry_type'), entry.get('id'), entry.get('accession'), self.description)



def family(id) -> PfamFamily:

    # Make the request
    root = request('/family/' + id)

    # Get the release version of the pfam database
    release_version = root.get('release')
    release_date = datetime.strptime(root.get('release_date'), '%Y-%m-%d')

    # Return the response parsed
    return PfamFamily(release_version, release_date, root[0])

def clan(id) -> PfamClan:

    # Make the request
    root = request('/clan/' + id)

    # Get the release version of the pfam database
    release_version = root.get('release')
    release_date = datetime.strptime(root.get('release_date'), '%Y-%m-%d')

    # Return the response parsed
    return PfamClan(release_version, release_date, root[0])

def families() -> [PfamEntry]:

    # Make the request
    text = request('/families', params={'output': 'text'})

    # Output list
    output: [PfamEntry] = []

    # Iterates over lines (family entry)
    for line in text.splitlines():

        # Split the line in three field
        fields = line.split('\t')
        
        if len(fields) == 3:
            output.append (PfamEntry('Pfam-A', fields[1], fields[0], fields[2]))

    # Return the response parsed
    return output

def clans() -> [PfamEntry]:

    # Make the request
    text = request('/clans', params={'output': 'text'})

    # Output list
    output: [PfamEntry] = []

    # Iterates over lines (family entry)
    for line in text.splitlines():

        # Split the line in three field
        fields = line.split('\t')
        
        if len(fields) == 3:
            output.append (PfamEntry('Clan', fields[1], fields[0], fields[2]))

    # Return the response parsed
    return output

def protein(id) -> PfamProtein:

    # Make the request
    root = request('/protein/' + id)

    # Get the release version of the pfam database
    release_version = root.get('release')
    release_date = datetime.strptime(root.get('release_date'), '%Y-%m-%d')

    # Return the response parsed
    return PfamProtein(release_version, release_date, root[0])


def proteins(family) -> [str]:

    # Options of the request
    params = {  'format': 'pfam',
                'alnType': 'full',
                'order': 'a',
                'case': 'l',
                'gaps': 'none',
                'download': '0',
                'output': 'text' }

    # Make the request
    text = request('/family/' + family + '/alignment/full/format', params=params)

    # Output list
    output: [str] = []

    # Iterates over lines
    for line in text.splitlines():

        # Split the line to get the name
        output.append(line.split('/')[0])

    # Return the list of protein names
    return output