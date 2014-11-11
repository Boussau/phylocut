//
// File: PhyloMerge.cpp
// Created by: Julien Dutheil, Bastien Boussau
// Created on: Sunday, December 2nd 2007 16:48
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>
#include <string.h>
#include <unordered_set>

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Distance.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/PhylipDistanceMatrixFormat.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Io/Newick.h>

const std::string ONE="E"; //Property used to store the taxon corresponding to parent and son0 nodes 
const std::string TWO="Fu"; //Property used to store the taxon corresponding to parent and son1 nodes 
const std::string THREE="S"; //Property used to store the taxon corresponding to son0 and son1 nodes 

//Compilation:
//g++  -pipe -o phylomerge bppPhyloSampler.cpp -I/usr/local/include  -L. -L/usr/local/lib  -g -Wall -fopenmp -std=c++0x -lbpp-core -lbpp-seq -lbpp-phyl


using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "phylomerge parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();

  (*ApplicationTools::message << "      Options considered: ").endLine();
(*ApplicationTools::message << " - input.method='tree' (could also be a distance matrix with the option 'matrix')").endLine();
(*ApplicationTools::message << " - input.tree.file=file with the tree in it.").endLine();
(*ApplicationTools::message << " - deletion.method='threshold' ou 'random' ou 'sample' ou 'taxon'. 'threshold' removes sequences that are so close in the tree that their distance is lower than the 'threshold' value (which is given as another option to the program, default is 0.01). 'sample': random choice of sample_size sequences (default is 10). 'taxon': choice is guided by the identity of the species the sequences come from. In cases several sequences from the same species are monophyletic, a choice will be made according to the 'choice.criterion' option").endLine();
(*ApplicationTools::message << " - choice.criterion='length' ou  'length.complete' ou 'merge'. 'length' means the longest sequence is selected. 'length.complete' : means the largest number of complete sites (no gaps). 'merge' means that the set of monophyletic sequences is used to build one long 'chimera' sequence corresponding to the merging of them.").endLine();
(*ApplicationTools::message << " - selection.by.taxon='no' ou 'yes'").endLine();
(*ApplicationTools::message << " - sequence.to.taxon=linkfile: format: sequence name: species name. Can be replaced by the option taxon.to.sequence").endLine();
(*ApplicationTools::message << " - taxon.to.sequence=linkfile: format: species name: sequence name. Can be replaced by the option sequence.to.taxon").endLine();
(*ApplicationTools::message << " - taxons.to.remove= file containing set of species from which sequences should be removed").endLine();
(*ApplicationTools::message << " - taxons.to.refine= file containing set of species on which the sampling/merging should be done. If not specified, all species are concerned.").endLine();
(*ApplicationTools::message << " - prescreening.on.size.by.taxon='no' : removes the sequences that are very short compared to other sequences of the same species. If there is only one sequence in this species, it is not removed.").endLine();
(*ApplicationTools::message << " - output.sequence.file=refined alignment").endLine();



  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}


/****************************************************************
 * Annotate nodes with taxon names when possible.
 ****************************************************************/
void postOrderAnnotateNodesWithTaxa(TreeTemplate< Node > &tree, Node *node, std::map<std::string, std::string>& seqSp)
{
  if (node->isLeaf()) 
  {
  //  node->setNodeProperty(ONE, BppString(seqSp[node->getName()]));
  //  node->setNodeProperty(TWO, BppString(seqSp[node->getName()]));
    
    if (tree.getRootNode() == node) //if at the root, we annotate the two properties including the (absent) father node
    {
      node->setNodeProperty(ONE, BppString(seqSp[node->getName()]));
      node->setNodeProperty(TWO, BppString(seqSp[node->getName()]));
      postOrderAnnotateNodesWithTaxa(tree, node->getSon(0), seqSp);
      node->setNodeProperty(THREE, *(node->getSon(0)->getNodeProperty(THREE) ) );
    }
    else {
      node->setNodeProperty(THREE, BppString(seqSp[node->getName()]));
    }
  }
  else 
  {
    vector<string> sonTaxa;
    for (unsigned int i = 0 ; i < node->getNumberOfSons() ; i++) 
    {
      postOrderAnnotateNodesWithTaxa(tree, node->getSon(i), seqSp);
      sonTaxa.push_back((dynamic_cast<const BppString*>(node->getSon(i)->getNodeProperty(THREE)))->toSTL());
    }
    vector<string> sonTaxaUnique = VectorTools::unique(sonTaxa);
    if (sonTaxaUnique.size() >1) 
    {
      node->setNodeProperty(THREE, BppString("#"));
    }
    else 
    {
      node->setNodeProperty(THREE, BppString( sonTaxaUnique[0] ));
    }
  }
  return;
}


/****************************************************************
 * Annotate nodes with taxon names when possible.
 ****************************************************************/
void preOrderAnnotateNodesWithTaxa(TreeTemplate< Node > &tree, Node *node, std::map<std::string, std::string> seqSp){
  if ((node->isLeaf()) && node->hasFather()) 
  { //A leaf that is not the root
    string fatherTaxa;
    //Here we assume bifurcation.
  if (node == node->getFather()->getSon(0))
  {
    fatherTaxa = ((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(TWO)))->toSTL());

  }
  else 
  {
    fatherTaxa = ((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(ONE)))->toSTL());
  }
    string sonTaxa = ((dynamic_cast<const BppString*>(node->getNodeProperty(THREE)))->toSTL());
    if (fatherTaxa == sonTaxa) 
    {
      node->setNodeProperty(ONE, BppString(sonTaxa));
      node->setNodeProperty(TWO, BppString(sonTaxa));
    }
    else 
    {
      node->setNodeProperty(ONE, BppString("#"));
      node->setNodeProperty(TWO, BppString("#"));
    }

  }
  if ((!node->isLeaf()) && node->hasFather()) 
  { //Not a leaf, not at the root
    string fatherTaxa;
    //Here we assume bifurcation.
    if (node == node->getFather()->getSon(0))
    {
/*      std::cout <<"node id: "<<node->getFather()->getId()<<std::endl;
      std::cout <<"root id: "<<tree.getRootNode()->getId()<<std::endl;
      std::cout <<"root is leaf?: "<<tree.getRootNode()->isLeaf()<<std::endl;
      std::cout <<"root numSons?: "<<tree.getRootNode()->getNumberOfSons()<<std::endl;*/
      fatherTaxa = ((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(TWO)))->toSTL());
    }
    else 
    {
      fatherTaxa = ((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(ONE)))->toSTL());
    }
    string son0Taxa = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(THREE)))->toSTL());
    string son1Taxa = ((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty(THREE)))->toSTL());
    if (fatherTaxa == son0Taxa) 
    {
      node->setNodeProperty(ONE, BppString(son0Taxa));
    }
    else 
    {
      node->setNodeProperty(ONE, BppString("#"));
    }
    if (fatherTaxa == son1Taxa) 
    {
      node->setNodeProperty(TWO, BppString(son1Taxa));
    }
    else 
    {
      node->setNodeProperty(TWO, BppString("#"));
    }
  }
  for (unsigned int i = 0 ; i < node->getNumberOfSons() ; i++) 
  {
    preOrderAnnotateNodesWithTaxa(tree, node->getSon(i), seqSp);
  }
  return;
}



/****************************************************************
 * Annotate all nodes of a tree with taxon names when possible.
 ****************************************************************/
void annotateTreeWithSpecies(TreeTemplate<Node>& tree, std::map<std::string, std::string>& seqSp) {
    
    //Now we have to go through the tree and remove sequences in groups of sequences from the same taxon
    //We do a double-recursive tree traversal to annotate nodes by taxa below them, in all three directions.
    //We use 3 node properties:
    //ONE: subtree defined by parent and son 0
    //TWO: subtree defined by parent and son 1
    //THREE: subtree defined by son 0 and son 1.
    
    postOrderAnnotateNodesWithTaxa(tree, tree.getRootNode(), seqSp);
    
    //Filling the 2 missing properties at the root
    Node* root = tree.getRootNode();
    if (!root->hasNodeProperty(TWO) )
    {
        
          /*       root->setNodeProperty(ONE, BppString("#"));
         root->setNodeProperty(TWO, BppString("#"));*/

        
        root->setNodeProperty(ONE, BppString(dynamic_cast<const BppString*>(root->getSon(0)->getNodeProperty(THREE) )->toSTL() ) );
        root->setNodeProperty(TWO, BppString(dynamic_cast<const BppString*>(root->getSon(1)->getNodeProperty(THREE) )->toSTL() ) );
        
        /*
         vector<string>     sonTaxa;
         //We assume bifurcation
         sonTaxa.push_back((dynamic_cast<const BppString*>(root->getSon(0)->getNodeProperty(THREE)))->toSTL());
         sonTaxa.push_back((dynamic_cast<const BppString*>(root->getSon(1)->getNodeProperty(THREE)))->toSTL());
         vector<string> sonTaxaUnique = VectorTools::unique(sonTaxa);
         if (sonTaxaUnique.size() >1)
         {
         root->setNodeProperty(ONE, BppString("#"));
         root->setNodeProperty(TWO, BppString("#"));
         }
         else
         {
         root->setNodeProperty(ONE, BppString( sonTaxaUnique[0] ));
         root->setNodeProperty(TWO, BppString( sonTaxaUnique[0] ));
         }*/
    }
    preOrderAnnotateNodesWithTaxa(tree, tree.getRootNode(), seqSp);
    //Now the tree has nodes annotated with taxon names
    Nhx *nhx = new Nhx();
    nhx->write(tree, cout);
    delete nhx;
}



double getBootstrapValueOnBranchBetweenTheseTwoNodes ( Node* it, Node* next ) {
    if (it->hasFather () && it->getFather () == next) {
        if (it->hasBootstrapValue()) {
            std::cout << "HAS BOOTSTRAP" <<std::endl;
        }
        return ( it->getBootstrapValue() );
    }
    else {
        if (next->hasBootstrapValue()) {
            std::cout << "HAS BOOTSTRAP" <<std::endl;
        }
        return ( next->getBootstrapValue() );
    }
    
}


/****************************************************************
 * Moves cutNode to be placed next to newBrother
 ****************************************************************/
void group(Node* newBrother, Node* cutNode, TreeTemplate<Node>& tree) {
    Node *oldFather, *oldGrandFather, *brother, *newBrothersFather, *N;
    double dist = 0.1;
    newBrothersFather = newBrother->getFather();
    oldFather = cutNode->getFather();
    
    
    //Get all old brothers ; a binary tree is assumed here (because of the "break")
    for(unsigned int i=0;i<oldFather->getNumberOfSons();i++)
        if(oldFather->getSon(i)!=cutNode){brother=oldFather->getSon(i); break;}
    
    if (!(oldFather->hasFather())) {//we displace the outgroup, need to reroot the tree
        //NB : brother is the other son of the root
        int id0 = oldFather->getId();
        int idBrother = brother->getId();
        N=new Node();
        
        N->addSon(newBrother);
        newBrother->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        
        //we remove cutNode from its old neighborhood
        for(unsigned int i=0;i<oldFather->getNumberOfSons();i++) {
            if(oldFather->getSon(i)==cutNode){oldFather->removeSon(i); break;}
        }
        // we move node cutNode
        N->addSon(cutNode);
        cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        
        // update N neighbours
        for(unsigned int i=0;i<newBrothersFather->getNumberOfSons();i++)
            if(newBrothersFather->getSon(i)==newBrother)
            {
                newBrothersFather->setSon(i, N); break;
            }
        N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        
        unsigned int oldFatherId = oldFather->getId();
        tree.rootAt(brother->getId());
        for(unsigned int i=0;i<brother->getNumberOfSons();i++) {
            if(brother->getSon(i)==oldFather){brother->removeSon(i);break;}
        }
        if(tree.hasNode(oldFatherId))
            delete oldFather;
        //We renumber the nodes
        brother->setId(id0);
        N->setId(idBrother);
    }
    else  {
        int id0 = oldFather->getId();
        //we create a new node N which will be the father of cutNode and newBrother
        N=new Node();
        
        N->addSon(newBrother);
        newBrother->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        // we move node cutNode
        N->addSon(cutNode);
        
        cutNode->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        
        // update N neighbours
        for(unsigned int i=0;i<newBrothersFather->getNumberOfSons();i++)
            if(newBrothersFather->getSon(i)==newBrother){newBrothersFather->setSon(i, N); break;}
        N->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        oldGrandFather = oldFather->getFather();
        for(unsigned int i=0;i<oldGrandFather->getNumberOfSons();i++)
            if(oldGrandFather->getSon(i)==oldFather){oldGrandFather->setSon(i, brother); break;}
        brother->setDistanceToFather(dist); // BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
        delete oldFather;
        N->setId(id0);
    }
    return ;

    
    
    
    
}






/****************************************************************
 * Simple function to tell whether a node belongs to a given species.
 ****************************************************************/
bool isNodeFromSpecies ( Node* node, string species) {
    string nodeTaxa1 = ((dynamic_cast<const BppString*>(node->getNodeProperty(ONE)))->toSTL());
    string nodeTaxa2 = ((dynamic_cast<const BppString*>(node->getNodeProperty(TWO)))->toSTL());
    string nodeTaxa3 = ((dynamic_cast<const BppString*>(node->getNodeProperty(THREE)))->toSTL());
    //std::cout << "species: "<<species << " nodeTaxa1 " << nodeTaxa1 << " nodeTaxa2 " <<nodeTaxa2 << " nodeTaxa3 " << nodeTaxa3  ;

    if ( nodeTaxa1 == species || nodeTaxa2 == species || nodeTaxa3 == species ) {
        return true;
    }
    else {
        return false;
    }
}


/****************************************************************
 * Finds sets of connected nodes that share the same origin species.
 ****************************************************************/

void searchForConnectedComponents(Node* node, vector< vector < Node* > >& connectedNodes, const string& species) ;

void findConnectedNodes(Node* node, vector< vector<Node*> >& connectedNodes, size_t i, const string species) {
    if (isNodeFromSpecies(node, species) ) {
        connectedNodes[i].push_back(node);
        if (node->getNumberOfSons() > 0 ) {
            vector<Node*> nodes = node->getSons();
            for (auto n = nodes.begin(); n != nodes.end(); ++n) {
                findConnectedNodes(*n, connectedNodes, i, species);
                
            }
        }
    }
    else {
        if (node->getNumberOfSons() > 0 ) {
            vector<Node*> nodes = node->getSons();
            for (auto n = nodes.begin(); n != nodes.end(); ++n) {
                searchForConnectedComponents(*n, connectedNodes, species);
                
            }
        }

    }
    
}


void searchForConnectedComponents(Node* node, vector< vector < Node* > >& connectedNodes, const string& species) {
    
    if (isNodeFromSpecies(node, species) ) {
        std::cout << "Found one" <<std::endl;
        vector<Node* > nodes;
        connectedNodes.push_back(nodes);
        findConnectedNodes( node, connectedNodes, connectedNodes.size() - 1, species);
    }
    else {
        if (node->getNumberOfSons() > 0 ) {
            vector<Node*> sons = node->getSons();
            for (auto n = sons.begin(); n != sons.end(); ++n) {
                searchForConnectedComponents(*n, connectedNodes, species);
            }
        }
    }
}

/****************************************************************
 * Reorganizes the tree to group sequences coming from the same species.
 ****************************************************************/
void groupSequencesFromSpecies( TreeTemplate< Node > &tree,
                               const string& species,
                               const double& bootstrapThreshold) {
    vector<Node*> nodes = tree.getNodes();
    vector<Node*> nodesOfSpecies ;
    for (auto node = nodes.begin(); node != nodes.end() ; ++node) {
        if (isNodeFromSpecies(*node, species) )
            nodesOfSpecies.push_back( *node );
    }
    //Now we have a subset of nodes from species "species".
    if (nodesOfSpecies.size() == 1) {
        return;
    }
    //Are they all connected?
    vector< vector< Node* > > connectedNodes;
    std::cout << "groupSequencesFromSpecies 1" <<std::endl;
    std::cout << "tree.getRootNode(): " << tree.getRootNode()->getId() << std::endl;
    searchForConnectedComponents(tree.getRootNode(), connectedNodes, species);
    std::cout << "groupSequencesFromSpecies 2" <<std::endl;

    //Now connectedNodes contains all connected components
    if (connectedNodes.size() == 1) {
        //They are all connected, we exit
        return;
    }
    //Several connected components, we need to see if we can group them.
    std::cout << "groupSequencesFromSpecies: connectedNodes.size(): "<< connectedNodes.size() <<std::endl;

    for (size_t i = 0; i < connectedNodes.size()-1 ; ++i) {
        for (size_t j = i; j < connectedNodes.size() ; ++j) {
            std::cout << "groupSequencesFromSpecies 21" <<std::endl;

            vector<Node*> path = TreeTemplateTools::getPathBetweenAnyTwoNodes(*(connectedNodes[i][0]), *(connectedNodes[j][0]) );
            std::cout << "groupSequencesFromSpecies 22" <<std::endl;

            //We have the path. We are only interested in nodes not annotated with the species of interest
            vector<Node*> pathBetweenComponents;
            size_t index = 0;
            for (auto it = path.begin(); it != path.end(); ++it) {
                if ( ! isNodeFromSpecies(*it, species) ) {
                    if (pathBetweenComponents.size() == 0) {
                        pathBetweenComponents.push_back(path[index-1]);
                    }
                    pathBetweenComponents.push_back(*it);
                }
                else if (pathBetweenComponents.size() > 0) {
                    pathBetweenComponents.push_back(path[index]);
                    break;
                }
                index++;
            }
            std::cout << "groupSequencesFromSpecies 3" <<std::endl;

            //Now we have the path between the components, with one node in the first component, and one in the last
            //Does any of these nodes include a high bootstrap branch, which would mean we can't join them?
            bool canJoinThem = true;
            for (auto it = pathBetweenComponents.begin(); it != pathBetweenComponents.end(); ++it) {
                if (it+1 == pathBetweenComponents.end() ) {
                    break;
                }
                Node* next = *(it+1);
                //Find the branch on which we want to check the bootstrap
                double bootstrap = getBootstrapValueOnBranchBetweenTheseTwoNodes ( *it, next );
                if (bootstrap > 1) {
                    bootstrap = bootstrap/100.0;
                }
                if (bootstrap > bootstrapThreshold) {
                    canJoinThem = false;
                    break;
                }
            }
            std::cout << "groupSequencesFromSpecies 4" <<std::endl;

            if (canJoinThem) {
                //Need to rearrange the tree.
                Node* node1 = *( pathBetweenComponents.begin() );
                Node* node2 = *( pathBetweenComponents.end()-1 );
                if ( node1->hasFather() && node1->getFather() == pathBetweenComponents[1]) {
                    group(node1, node2, tree);
                }
                else {
                    group(node2, node1, tree);
                }
            }
            std::cout << "groupSequencesFromSpecies 5" <<std::endl;

        }
    }
    
}


/****************************************************************
 * Refines the tree to group sequences according to their species in situations where
 * there is no highly supported edge that separates the said sequences.
 ****************************************************************/
void refineTree ( TreeTemplate< Node > &tree,
                 const std::map<std::string, std::string>& seqSp,
                 const unordered_set<string>& speciesToRefine,
                 const double& bootstrapThreshold) {
    std::cout << "Refining the tree..." <<std::endl;
    if (speciesToRefine.size() != 0 ) {//Then all species need working with
        std::cout << "in refineTree 2" <<std::endl;
        for (auto it = speciesToRefine.begin(); it != speciesToRefine.end(); ++it) {
            std::cout << "Species: "<< *it <<std::endl;
            if (*it != "")
                groupSequencesFromSpecies (tree, *it, bootstrapThreshold);
        }
        std::cout << "in refineTree 21" <<std::endl;

        
    }
    else {
        std::set<string> species;
        std::cout << "in refineTree 3" <<std::endl;

        for (auto it = seqSp.begin(); it!= seqSp.end(); ++it) {
            species.insert(it->second);
        }
        std::cout << "in refineTree 4" <<std::endl;

        for (auto it = species.begin(); it != species.end(); ++it) {
            std::cout << "Species: "<< *it <<std::endl;
            if (*it != "")
                groupSequencesFromSpecies (tree, *it, bootstrapThreshold);
        }
        std::cout << "in refineTree 41" <<std::endl;

        
    }
    
}



/****************************************************************
 * Merges sequences to produce one new sequence.
 ****************************************************************/
string buildMergedSequence(vector<string>& descendantSequences, const VectorSiteContainer& seqs) {
	vector<string> uniqueSequences = VectorTools::unique(descendantSequences);
	/*VectorTools::print (descendantSequences);
		VectorTools::print (uniqueSequences);*/
	/*
	string sequence = seqs.getSequence(descendantSequences[0]).toString();
    for (unsigned int i = 0; i < seqs.getNumberOfSites(); i++)
    {
        const int* element = &seqs(seqs.getSequencePosition (descendantSequences[0]), i);
        if (seqs.getAlphabet()->isGap(*element) || seqs.getAlphabet()->isUnresolved(*element) ) {
            for (unsigned int j = 1; j < descendantSequences.size() ; j++)
            {
                const int* element2 = &seqs(seqs.getSequencePosition (descendantSequences[j]), i);
                if ( ! ( seqs.getAlphabet()->isGap(*element2) || seqs.getAlphabet()->isUnresolved(*element2) ) )
                    //Put the site in the right place in the string sequence
                    sequence.replace(i, 1, seqs.getAlphabet()->intToChar(*element2));
                   // sequence[i] = string( ( seqs.getAlphabet()->intToChar(*element2) ) ); 
                    break;
                    }
        }
    }  */
  VectorSequenceContainer* selSeqs = new VectorSequenceContainer(seqs.getAlphabet());
	//SequenceContainer* selSeqs = 0;
	SequenceContainerTools::getSelectedSequences (seqs, descendantSequences, *selSeqs, true);
	//use consensus instead:
	VectorSiteContainer *selSeqsSites = new VectorSiteContainer ( *selSeqs );
	//		Sequence* sequence = SiteContainerTools::getConsensus(*(dynamic_cast <const SiteContainer*> (selSeqs) ), "consensus", true, false);
	Sequence* sequence = SiteContainerTools::getConsensus(*selSeqsSites, "consensus", true, false);
	string toReturn = sequence->toString();
	delete selSeqs;
	delete selSeqsSites;
	delete sequence;
    return toReturn ;
}


/****************************************************************
 * Select only one sequence in a vector of sequence names.
 ****************************************************************/
string selectSequenceAmongSequences(vector<string>& descendantSequences, string critMeth, VectorSiteContainer& seqs, string name = "") 
{
  if (critMeth == "length" || critMeth == "length.complete")
  {
    vector<double> seqLen;
    for(unsigned int i = 0; i < descendantSequences.size(); i++)
    {
      if(critMeth == "length.complete")
        seqLen.push_back(SequenceTools::getNumberOfCompleteSites(seqs.getSequence(descendantSequences[i])));
      else
        seqLen.push_back(SequenceTools::getNumberOfSites(seqs.getSequence(descendantSequences[i])));
    }      
    return descendantSequences[VectorTools::whichMax(seqLen)];
  }
  else if (critMeth == "random")
  {
    return descendantSequences[RandomTools::giveIntRandomNumberBetweenZeroAndEntry (descendantSequences.size())];
  }
  else if (critMeth == "merge") {
      string seq = buildMergedSequence(descendantSequences, seqs);
      BasicSequence bseq = BasicSequence (name+"_"+descendantSequences[0], seq, seqs.getAlphabet());
      //MYSTERY
       if (!seqs.hasSequence(name+"_"+descendantSequences[0])) seqs.addSequence(bseq);
      return name+"_"+descendantSequences[0];
  }
  else throw Exception("Unknown criterion: " + critMeth);  
}


/****************************************************************
 * Build a map containing distances between ancestor node and leaves.
 ****************************************************************/
map <string, double> computeDistanceBetweenLeavesAndAncestor(TreeTemplate< Node > &tree, Node &node, vector <string>& sameTaxaSequences) {
    map <string, double> sequenceDistances;
    for (unsigned int i = 0 ; i < sameTaxaSequences.size(); i ++) {
        sequenceDistances.insert( pair<string,double>(sameTaxaSequences[i], TreeTemplateTools::getDistanceBetweenAnyTwoNodes (node, *(tree.getNode(sameTaxaSequences[i]) ) ) ) );
    }
    return sequenceDistances;
}


/****************************************************************
 * Returns a vector with strings corresponding to entries in a 
 * string/double map whose values match a particular value.
 ****************************************************************/
vector<string> getSequencesAtAGivenDistanceFromAncestor( map <string, double>& sequenceDistances, double value) {
    map<string,double>::iterator it;
    vector<string> sequences;
    for ( it = sequenceDistances.begin() ; it != sequenceDistances.end(); it++ ) {
        if( (*it).second == value ) {
            sequences.push_back( (*it).first );
        }
    }
    return sequences;
}


/****************************************************************
 * Sort a string vector given a map linking the strings to doubles.
 ****************************************************************/
void sortVectorGivenMap(vector <string>& sameTaxaSequences, map <string, double>& sequenceDistances){
    vector <double> values;
    sameTaxaSequences.clear();
    map<string,double>::iterator it;
    for ( it = sequenceDistances.begin() ; it != sequenceDistances.end(); it++ ) {
        values.push_back( (*it).second );   
    }
    sort (values.begin(), values.end() );
    for (unsigned int i = 0 ; i < values.size() ; i++) {
        VectorTools::append( sameTaxaSequences, getSequencesAtAGivenDistanceFromAncestor(sequenceDistances , values[i]) );
    }
    return;
}


/****************************************************************
 * Tree Traversal to select sequences non-redundant in terms of taxa: 
 * when there is monophyly, we select/build only one sequence.
 ****************************************************************/
vector<string> selectSequencesToKeep(TreeTemplate< Node > &tree, 
                                     Node *node, 
                                     vector <string>& sequencesAncestor, 
                                     string critMeth, 
                                     VectorSiteContainer& seqs, 
                                     unordered_set<string> speciesToRefine )/*,
                                     const map <string, double>& sequenceDistances)*/
{

//std::cout << TreeTemplateTools::treeToParenthesis( tree, true  ) <<std::endl;

  vector<string> selectedSeqs;
  string nodeTaxa1 = ((dynamic_cast<const BppString*>(node->getNodeProperty(ONE)))->toSTL()); 
  string nodeTaxa2 = ((dynamic_cast<const BppString*>(node->getNodeProperty(TWO)))->toSTL()); 
  string nodeTaxa3 = ((dynamic_cast<const BppString*>(node->getNodeProperty(THREE)))->toSTL()); 
  vector<string> seqs0;
  vector<string> seqs1;
  map <string, double> sequenceDistances;

    //  std::cout << "BEFORE all: "<< nodeTaxa1 << " and " << nodeTaxa2 << " and " << nodeTaxa3 << std::endl;


  if (! node->isLeaf()) 
  { //If not at a leaf
    //Here we assume bifurcation.
    seqs0 = TreeTemplateTools::getLeavesNames(*(node->getSon(0)));
    seqs1 = TreeTemplateTools::getLeavesNames(*(node->getSon(1)));
  }
  else {
  }
  if (node->hasFather()) 
  { //Not at the root
    if (nodeTaxa1 != "#") 
    { //If father node and son0 node come from the same taxon
      bool sampleHere = true;
      if (! node->isLeaf()) 
      {
        string nodeSon0Taxa1 = ((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty(ONE)))->toSTL());
        string nodeSon0Taxa2 = ((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty(TWO)))->toSTL());
        if (nodeSon0Taxa1 == nodeTaxa1 || nodeSon0Taxa2 == nodeTaxa1) {
          sampleHere = false;//We can sample at son1 node instead
        }
      }
      else 
      {
        // selectedSeqs.push_back(node->getName());
        // sequencesAncestor.push_back(node->getName());
          
       //TEMP sampleHere = false;
      }
      if (sampleHere) 
      {       //we sample here
        vector <string> sameTaxaSequences = VectorTools::vectorUnion(sequencesAncestor, seqs0);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, sameTaxaSequences);
        //we sort the vector of sequences with respect to the distance to the ancestor
        sortVectorGivenMap(sameTaxaSequences, sequenceDistances);
	if (speciesToRefine.count(nodeTaxa1) != 0 || speciesToRefine.size()==0 ) {
        	if (critMeth == "merge")
          		selectedSeqs.push_back(selectSequenceAmongSequences(sameTaxaSequences, critMeth, seqs, nodeTaxa1));
        	else 
          		selectedSeqs.push_back(selectSequenceAmongSequences(sameTaxaSequences, critMeth, seqs));
	}
	else {
		VectorTools::append(selectedSeqs, sameTaxaSequences);
	}
          //We have to go further in the unsampled tree
          vector <string> temp = VectorTools::vectorUnion(sequencesAncestor, seqs0);
          sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
          sortVectorGivenMap(temp, sequenceDistances);
          VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), temp, critMeth, seqs, speciesToRefine) );
        
      }
      else if (! node->isLeaf() ) {
        //Or we do not sample here, and we have to go further in the tree on both sides
        vector <string> temp = VectorTools::vectorUnion(sequencesAncestor, seqs0);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), temp, critMeth, seqs, speciesToRefine) );
        temp = VectorTools::vectorUnion(sequencesAncestor, seqs1);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), temp, critMeth, seqs, speciesToRefine) );
        
      }
    }
    
    if (nodeTaxa2 != "#") 
    {      
      bool sampleHere = true;
      if (! node->isLeaf()) 
      {
        string nodeSon0Taxa1 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(ONE)))->toSTL());
        string nodeSon0Taxa2 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(TWO)))->toSTL());
        if (nodeSon0Taxa1 == nodeTaxa2 || nodeSon0Taxa2 == nodeTaxa2) {
          sampleHere = false; //We can sample at son0 node instead
        }
      }
      else 
      {
        // selectedSeqs.push_back(node->getName());
        // sequencesAncestor.push_back(node->getName());
       //TEMP sampleHere = false;
      }
      if (sampleHere) 
      { //we sample here
        vector <string> sameTaxaSequences = VectorTools::vectorUnion(sequencesAncestor, seqs1);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, sameTaxaSequences);
        sortVectorGivenMap(sameTaxaSequences, sequenceDistances);

	if (speciesToRefine.count(nodeTaxa2) != 0 || speciesToRefine.size()==0 ) {
        	if (critMeth == "merge")
          		selectedSeqs.push_back(selectSequenceAmongSequences(sameTaxaSequences, critMeth, seqs, nodeTaxa2));
        	else 
          		selectedSeqs.push_back(selectSequenceAmongSequences(sameTaxaSequences, critMeth, seqs));
	}
	else {
		VectorTools::append(selectedSeqs, sameTaxaSequences);
	}
        //We have to go further in the unsampled tree
        vector <string> temp = VectorTools::vectorUnion(sequencesAncestor, seqs1);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), temp, critMeth, seqs, speciesToRefine) );
      }
      else if (! node->isLeaf() ) {
        //Or we do not sample here, and we have to go further in the tree on both sides
        vector <string> temp = VectorTools::vectorUnion(sequencesAncestor, seqs0);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), temp, critMeth, seqs, speciesToRefine) );
        temp = VectorTools::vectorUnion(sequencesAncestor, seqs1);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), temp, critMeth, seqs, speciesToRefine) );
        
      }
    }
    
    if (nodeTaxa3 != "#") 
    {      
      bool sampleHere = true;
      if (node->hasFather()) 
      {
        string fatherTaxa3 = ((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(THREE)))->toSTL());
        if (fatherTaxa3 !="#" ) 
        {
          sampleHere = false; //We can sample at the father node instead (hopefully, we have)
        }
      }
      if (sampleHere ) 
      {
        vector <string> descendantSequences;
        if (! node->isLeaf()) 
        {
          descendantSequences = VectorTools::vectorUnion(seqs0, seqs1);
        }
        else 
        {
          if (node->hasFather() ) 
          {
            bool sample = true;
            if (node->getFather()->getSon(0) == node) 
            {
              if (((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(ONE)))->toSTL()) != "#" && node->getFather()->hasFather()) 
                sample=false;
            }
            else 
            {
              if (((dynamic_cast<const BppString*>(node->getFather()->getNodeProperty(TWO)))->toSTL()) != "#"  && node->getFather()->hasFather()) 
                sample=false;
            }
            if (sample) 
            {
              descendantSequences.push_back(node->getName());
            }
          }
        }
        if (descendantSequences.size() > 0) {
          if (descendantSequences.size() > 1) {
            sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, descendantSequences);
            sortVectorGivenMap(descendantSequences, sequenceDistances);
          }

	if (speciesToRefine.count(nodeTaxa3) != 0 || speciesToRefine.size()==0 ) {
        	if (critMeth == "merge")
          		selectedSeqs.push_back(selectSequenceAmongSequences(descendantSequences, critMeth, seqs, nodeTaxa3));
        	else 
          		selectedSeqs.push_back(selectSequenceAmongSequences(descendantSequences, critMeth, seqs));
	}
	else {
		VectorTools::append(selectedSeqs, descendantSequences);
	}

          
        }
      }
    }
    else if (nodeTaxa1 == "#" && nodeTaxa2 == "#" && nodeTaxa3 == "#") {
      if (! node->isLeaf() ) {
        //We do not sample here, and we have to go further in the tree on both sides        
        vector <string> temp = VectorTools::vectorUnion(sequencesAncestor, seqs0);                
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), temp, critMeth, seqs, speciesToRefine) );
        temp = VectorTools::vectorUnion(sequencesAncestor, seqs1);
        sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
        sortVectorGivenMap(temp, sequenceDistances);
        VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), temp, critMeth, seqs, speciesToRefine) );
        
      }
    }
    else {/*
    Nhx *nhx = new Nhx();
    nhx->write(tree, cout);    

      if (! node->isLeaf() ) 
        std::cerr<< "Problem: a node is not properly annotated with taxon names or #. This is node: "<< node->getId() << " in tree " << TreeTemplateTools::treeToParenthesis( tree, true  ) <<std::endl;
      else 
        std::cerr<< "Problem: a node is not properly annotated with taxon names or #. This is node: "<< node->getName() << " in tree " << TreeTemplateTools::treeToParenthesis( tree, true  ) <<std::endl;
*/
    }
  }
  else
  { //At the root
    //  std::cout << "At ThE root: "<< nodeTaxa1 << " and " << nodeTaxa2 << " and " << nodeTaxa3 << std::endl;
    if (nodeTaxa1 == nodeTaxa2 && nodeTaxa1 == nodeTaxa3 && nodeTaxa1 != "#") 
    { //All the sequences are from the same taxon
      vector <string> temp = VectorTools::vectorUnion(seqs0, seqs1); 
      temp.push_back(node->getName());
      sequenceDistances = computeDistanceBetweenLeavesAndAncestor(tree, *node, temp);
      sortVectorGivenMap(temp, sequenceDistances);
	if (speciesToRefine.count(nodeTaxa3) != 0 || speciesToRefine.size()==0 ) {
        	if (critMeth == "merge")
          		selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeTaxa3));
        	else 
          		selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
	}
	else {
            VectorTools::append(selectedSeqs, temp);
        }
    }
    else //Not all the sequences are from the same taxon
    {
      //we may have to sample the root sequence
        if (node->isLeaf()) //Root node is a leaf
        {
        string nodeSon0Taxa1 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(ONE)))->toSTL());
        string nodeSon0Taxa2 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(TWO)))->toSTL());
            if ( nodeTaxa1 != nodeSon0Taxa1 && nodeTaxa1 != nodeSon0Taxa2 ) 
            { //We have to sample the root node
            vector <string> temp ;
            temp.push_back(node->getName());

        if (speciesToRefine.count(nodeTaxa1) != 0 || speciesToRefine.size()==0 ) {
                if (critMeth == "merge")
                    selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeTaxa1));
                else 
                    selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
        }
        else {
            VectorTools::append(selectedSeqs, temp);
        }
            }
        sequencesAncestor.push_back(node->getName());
                
        if (node->isLeaf()) {
            VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));
        }
        else {
            //we assume bifurcations
            VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));
            VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), sequencesAncestor, critMeth, seqs, speciesToRefine));
            }
        }
      else {//Root Node is not a leaf
        string nodeSon0Taxa1 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(ONE)))->toSTL());
        string nodeSon0Taxa2 = ((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty(TWO)))->toSTL());
        string nodeSon1Taxa1 = ((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty(ONE)))->toSTL());
        string nodeSon1Taxa2 = ((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty(TWO)))->toSTL());
      //  std::cout << "nodeSon0Taxa1: "<< nodeSon0Taxa1 << " and " << nodeSon0Taxa2 << " and " << nodeSon1Taxa1 << " and " << nodeSon1Taxa2 <<std::endl;
        
        if (nodeSon0Taxa1 != "#" || nodeSon0Taxa2 != "#" ||nodeSon1Taxa1 != "#" || nodeSon1Taxa2 != "#" )
        {
         //We should group together some sequences on each side of the root 
            if (nodeSon0Taxa1 != "#") {
                vector <string> temp ;
                if (node->getSon(0)->isLeaf()){
                    temp.push_back(node->getSon(0)->getName());
                }
                else {
                    VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(0)->getSon(0)))  );    
                }
                VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(1)))  );
                VectorTools::print(temp);
                if (speciesToRefine.count(nodeSon0Taxa1) != 0 || speciesToRefine.size()==0 ) {
                    if (critMeth == "merge")
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeSon0Taxa1));
                    else 
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
                }
                else {
                    VectorTools::append(selectedSeqs, temp);
                }
                //Continue the recursion
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0)->getSon(1), sequencesAncestor, critMeth, seqs, speciesToRefine));

            }
            else if (nodeSon0Taxa2 != "#") {
                vector <string> temp ;
                if (node->getSon(0)->isLeaf()){
                    temp.push_back(node->getSon(0)->getName());
                }
                else {
                    VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(0)->getSon(1)) ) );
                }
                VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(1)))  );

                                VectorTools::print(temp);

                if (speciesToRefine.count(nodeSon0Taxa2) != 0 || speciesToRefine.size()==0 ) {
                    if (critMeth == "merge")
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeSon0Taxa2));
                    else 
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
                }
                else {
                    VectorTools::append(selectedSeqs, temp);
                }
                //Continue the recursion
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0)->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));
            }
            else if (nodeSon1Taxa1 != "#") {
                vector <string> temp ;
                if (node->getSon(1)->isLeaf()){
                    temp.push_back(node->getSon(1)->getName());
                }
                else {
                    VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(1)->getSon(0))) );
                }
                VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(0)))  );

                                VectorTools::print(temp);

                if (speciesToRefine.count(nodeSon1Taxa1) != 0 || speciesToRefine.size()==0 ) {
                    if (critMeth == "merge")
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeSon1Taxa1));
                    else 
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
                }
                else {
                    VectorTools::append(selectedSeqs, temp);
                }
                //Continue the recursion
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1)->getSon(1), sequencesAncestor, critMeth, seqs, speciesToRefine));
            }
            else if (nodeSon1Taxa2 != "#") {
                vector <string> temp ;
                if (node->getSon(1)->isLeaf()){
                    temp.push_back(node->getSon(1)->getName());
                }
                else {
                    VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(1)->getSon(1))) );
                }                
                VectorTools::append(temp, TreeTemplateTools::getLeavesNames( *(node->getSon(0)))  );
                                VectorTools::print(temp);

                if (speciesToRefine.count(nodeSon1Taxa2) != 0 || speciesToRefine.size()==0 ) {
                    if (critMeth == "merge")
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeSon1Taxa2));
                    else 
                        selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
                }
                else {
                    VectorTools::append(selectedSeqs, temp);
                }
                //Continue the recursion
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1)->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));

            }
        }
        else {
            if (node->isLeaf()) 
            {
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));
            }
            else {
            //we assume bifurcations
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(0), sequencesAncestor, critMeth, seqs, speciesToRefine));
                VectorTools::append(selectedSeqs, selectSequencesToKeep(tree, node->getSon(1), sequencesAncestor, critMeth, seqs, speciesToRefine));
            }
        }
        
        /*if ( nodeTaxa1 != nodeSon0Taxa1 && nodeTaxa1 != nodeSon0Taxa2 && nodeTaxa1 != nodeSon1Taxa1 && nodeTaxa1 != nodeSon1Taxa2 && nodeTaxa1 != "#") 
        { //We have to sample the root node
          vector <string> temp ;
          temp.push_back(node->getName());
            if (speciesToRefine.count(nodeTaxa1) != 0 || speciesToRefine.size()==0 ) {
                if (critMeth == "merge")
                    selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs, nodeTaxa1));
                else 
                    selectedSeqs.push_back(selectSequenceAmongSequences(temp, critMeth, seqs));
            }
            else {
                VectorTools::append(selectedSeqs, temp);
            }
        }*/
      }


    }
  }
  return selectedSeqs;
}




class Index {
  public:
    double distance;
    unsigned int i1, i2;

  public:
    Index(double dist, unsigned int i, unsigned int j) : distance(dist), i1(i), i2(j) {}

  public:
    bool operator==(const Index& index) const { return distance == index.distance; }
    bool operator<(const Index& index) const { return distance < index.distance; }
};

class Test {
  private:
    unsigned int pos_;

  public:
    Test(unsigned int pos) : pos_(pos) {}
    
  public:
    bool operator()(const Index& index) { return index.i1 == pos_ || index.i2 == pos_; }
};


string replaceChar(string &str, const char ch1, const char ch2) {
  for (unsigned int i = 0; i < str.length(); ++i) {
    if (str[i] == ch1)
      str[i] = ch2;
  }

  return str;
}





int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*        Bio++ Phylogenetic Merger/Sampler, version 0.2          *" << endl;
  cout << "* Author: J. Dutheil, B. Boussau            Last Modif. 21/02/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication phylomerge(args, argv, "PhyloMerge");
  phylomerge.startTimer();

  //Get sequences:
  Alphabet* alphabet      = SequenceApplicationTools::getAlphabet(phylomerge.getParams());
  VectorSiteContainer* seqs = SequenceApplicationTools::getSiteContainer(alphabet, phylomerge.getParams());
  ApplicationTools::displayResult("Initial number of sequences:", seqs->getNumberOfSequences());

    
  string inputMethod = ApplicationTools::getStringParameter("input.method", phylomerge.getParams(), "tree");
  ApplicationTools::displayResult("Input method", inputMethod);

  DistanceMatrix* dist = 0;
  //  DistanceMatrix* matrix = 0;
  TreeTemplate< Node > * tree =0;
  if(inputMethod == "tree")
  {
      string treePath = ApplicationTools::getAFilePath("input.tree.file", phylomerge.getParams(), true, true);

      
      ifstream file_stream (treePath.c_str());
      vector<string> trees;
      int tree_i=0;
      if (file_stream.is_open())  //  ########## read trees ############
      {
          while (! file_stream.eof())
          {
              string line;
              getline (file_stream,line);
              if (line.find("(")!=line.npos)
              {
                  tree_i++;
                  trees.push_back(line);
              }
          }
      }
      if (tree_i > 1) {
          std::cout << "More than 1 tree in the input tree file, only considering the first tree." <<std::endl;
      }
      
      /*
	  TreeTemplate<Node> *treeTemp = dynamic_cast<TreeTemplate<Node>*> (PhylogeneticsApplicationTools::getTree(phylomerge.getParams() ) );
      if (treeTemp->getRootNode()->hasBootstrapValue()) {
          std::cout << "HAS BOOTSTRAP" <<std::endl;
      }
      else       {
          std::cout << "treetemp DOES NOT HAVE BOOTSTRAP" <<std::endl;
      }
tree = new TreeTemplate<Node>(*(treeTemp));
       delete treeTemp;*/
      tree = TreeTemplateTools::parenthesisToTree(trees[0], true);
      std::cout <<TreeTemplateTools::treeToParenthesis( *tree, true )<<std::endl;
      if (tree->getRootNode()->hasBootstrapValue()) {
          std::cout << "tree HAS BOOTSTRAP" <<std::endl;
      }
      else       {
          std::cout << "tree DOES NOT HAVE BOOTSTRAP" <<std::endl;
      }
    dist = TreeTemplateTools::getDistanceMatrix(*tree);
    
    //    VectorTools::print(tree->getLeavesNames());
  }
  else if(inputMethod == "matrix")
  {
    string distPath = ApplicationTools::getAFilePath("input.matrix", phylomerge.getParams(), true, true);
    PhylipDistanceMatrixFormat matIO;
    dist = matIO.read(distPath);
  }
  else throw Exception("Unknown input method: " + inputMethod);

  string deleteMeth = ApplicationTools::getStringParameter("deletion.method", phylomerge.getParams(), "threshold");
  ApplicationTools::displayResult("Deletion method", deleteMeth);

  string critMeth = ApplicationTools::getStringParameter("choice.criterion", phylomerge.getParams(), "length");
  ApplicationTools::displayResult("Sequence choice criterion", critMeth);

  bool useSpecies = ApplicationTools::getBooleanParameter("selection.by.taxon", phylomerge.getParams(), false, "", true, false);
  string speciesFile ;
  //We use a std::map to record the links between sequence name and species name.
  std::map<std::string, std::string> seqSp;

  unordered_set<string> speciesToRefine;

  if (useSpecies) 
  {
    if (ApplicationTools::getAFilePath("sequence.to.taxon", phylomerge.getParams(), true, true) != "none" ) {
    speciesFile = ApplicationTools::getAFilePath("sequence.to.taxon", phylomerge.getParams(), true, true);
    //In this file, the format is expected to be as follows :
    /*
     sequence1:SpeciesA
     sequence5:SpeciesA
     sequence3:SpeciesB
     ...
     */
    std::ifstream inSpSeq (speciesFile.c_str());
    std::string line;
    while(getline(inSpSeq,line)) {
		if (TextTools::hasSubstring (line, ":")) {
			//We divide the line in 2 : first, the sequence name, second the species name
			StringTokenizer st1 (line, ":", true);
			//Then we divide the sequence names
			seqSp.insert(make_pair(st1.getToken(0),st1.getToken(1)));
		}
    }
    deleteMeth="taxon";
    }
else {
    speciesFile = ApplicationTools::getAFilePath("taxon.to.sequence", phylomerge.getParams(), true, true);
    //In this file, the format is expected to be as follows :
    /*
     SpeciesA:sequence1
     SpeciesA:sequence5
     SpeciesB:sequence3
     ...
     */
    std::ifstream inSpSeq (speciesFile.c_str());
    std::string line;
    while(getline(inSpSeq,line)) {
		if (TextTools::hasSubstring (line, ":")) {
			//We divide the line in 2 : first, the species name, then the sequence name
			StringTokenizer st1 (line, ":", true);
			//Then we divide the sequence names
			seqSp.insert(make_pair(st1.getToken(1),st1.getToken(0)));
		}
    }
    deleteMeth="taxon";
}


    ///////////////////////////////////////////////////  
    //removing sequences from species we do not want
    ///////////////////////////////////////////////////
    string speciesToRemoveFile = ApplicationTools::getAFilePath("taxons.to.remove", phylomerge.getParams(), false, true);
    vector<string> speciesToRemove;
    std::string line;
    if (speciesToRemoveFile != "none") {
      std::ifstream inSp (speciesToRemoveFile.c_str());
      while(getline(inSp,line)) {
	line=TextTools::removeLastWhiteSpaces(line);
	line=TextTools::removeLastNewLines(line);
	replaceChar(line, ' ', '_');
	speciesToRemove.push_back(line);
      }
    }

    map<string,string>::iterator it;
    vector <string> sequencesToRemove;
      for ( it = seqSp.begin() ; it != seqSp.end(); it++ ) {
          //bool found = false;
          for (unsigned int i = 0 ; i<speciesToRemove.size() ; i++) {
              if (TextTools::hasSubstring((*it).first, speciesToRemove[i])) {
                  sequencesToRemove.push_back( (*it).first );
                  //found = true;
                  break;
              }
          }
          
          /*      if (!found) {
           std::cout << "not found "<< (*it).first <<std::endl;                                                                                                                                                    
           }*/
      }


    //Now we have sequences to remove
    vector <string> seqNames = seqs->getSequencesNames();

    vector <string> seqsToKeep ;

    vector <string> seqsToRemove;
	  
      if (sequencesToRemove.size() > 0 ) {
          VectorTools::diff(seqNames, sequencesToRemove, seqsToKeep);
          
          for (unsigned int i = 0 ; i < sequencesToRemove.size() ; i++) {
              std::cout <<"Removing sequence "<<sequencesToRemove[i] <<" # "<<i<<std::endl;
              //TreeTemplateTools::dropLeaf(*tree, sequencesToRemove[i]);
              if (seqs->hasSequence(sequencesToRemove[i]) ) {
                  seqs->deleteSequence(sequencesToRemove[i]);
                  seqSp.erase(sequencesToRemove[i]);
              }
              
          }
      }
    string nomi, nomj;
      
      //Instead we prune leaves
      for (unsigned int i = 0 ; i < sequencesToRemove.size() ; i++) {
          TreeTemplateTools::dropLeaf(*tree, sequencesToRemove[i]);
      }
      sequencesToRemove.clear();



    /////////////////////////////////////////////////// 
    //Selecting sequences from which we want to do the selection
    /////////////////////////////////////////////////// 
    string speciesToRefineFile = ApplicationTools::getAFilePath("taxons.to.refine", phylomerge.getParams(), false, true);
    if (speciesToRefineFile != "none") {
      std::ifstream inSp (speciesToRefineFile.c_str());
      while(getline(inSp,line)) {
	line=TextTools::removeLastWhiteSpaces(line);
	line=TextTools::removeLastNewLines(line);
	replaceChar(line, ' ', '_');
	speciesToRefine.emplace(line);
      }
    }
    // Now we have a list of species on which we will apply our algorithms



    ///////////////////////////////////////////////////
    //This pre-screening step aims at removing sequences that are short compared to other sequences for a given taxon
    //If there is only one sequence for a taxon, this sequence is not removed.
    ///////////////////////////////////////////////////
    bool preScreening = ApplicationTools::getBooleanParameter("prescreening.on.size.by.taxon", phylomerge.getParams(), false, "", true, false);
    if (preScreening) {
      ApplicationTools::displayTask("Pre-screening phase", true);  
      //First we build a map from species names to sequence names.
      //For one species name, we can have several sequence names
      std::map<std::string, std::deque<std::string> > spSeq;
      //      map<string,string>::iterator it;
      map<string, std::deque<std::string> >::iterator it2;
      for ( it = seqSp.begin() ; it != seqSp.end(); it++ ) {
        it2 = spSeq.find((*it).second);
        if ( it2 != spSeq.end() ) {
          (*it2).second.push_back((*it).first);
        }
        else {
          std::deque<string> temp (1, (*it).first);
          spSeq.insert( make_pair( (*it).second, temp) );
        }
      }
      //Now spSeq is filled
      //We go through each species. If we find sequences much smaller (1/2) than the largest one, we discard it
      vector <string> toKeep ; 
      for ( it2 = spSeq.begin() ; it2 != spSeq.end(); it2++ ) {
        vector<unsigned int> sizes;
        for (unsigned int i = 0 ; i < (*it2).second.size() ; i++) {
          BasicSequence seq = seqs->getSequence((*it2).second[i]);
          int tot = SequenceTools::getNumberOfSites(seq);
          sizes.push_back(tot);
        }
        unsigned int max = VectorTools::max(sizes);
        for (unsigned int i = 0 ; i < (*it2).second.size() ; i++) {
         // if (sizes[i] > max / 2 ) {
          if (sizes[i] > 50 && sizes[i] > max / 2) {
            toKeep.push_back((*it2).second[i]);
          }
          else {
              sequencesToRemove.push_back((*it2).second[i]);
              ApplicationTools::displayWarning("Discarding sequence " + (*it2).second[i] );
          }
        }
      }
      //Now we remove the sequences, if some are to remove...
      if (toKeep.size() != seqs->getNumberOfSequences()) {
          //Instead we prune leaves
          for (unsigned int i = 0 ; i < sequencesToRemove.size() ; i++) {
              TreeTemplateTools::dropLeaf(*tree, sequencesToRemove[i]);
          }
        //and we remove from the sequences        
        VectorSiteContainer asc(alphabet);
	  //SOME MYSTERY TO INVESTIGATE WITH FAMILY102
        for(unsigned int i = 0; i < toKeep.size(); i++)
	  if (!asc.hasSequence(toKeep[i])) asc.addSequence(seqs->getSequence(toKeep[i]));//MYSTERY
	//        asc.addSequence(seqs->getSequence(toKeep[i]));
        delete seqs;
        seqs = dynamic_cast<VectorSiteContainer* > ( asc.clone() );
      }
    }
  }

  string name;
  //Compute lengths:
  vector<string> seqNames;
  vector<unsigned int> seqLen(dist->size());
  for(unsigned int i = 0; i < dist->size(); i++)
  {
    name = dist->getName(i);
    if(critMeth == "length.complete")
      seqLen[i] = SequenceTools::getNumberOfCompleteSites(seqs->getSequence(name));
    else
      seqLen[i] = SequenceTools::getNumberOfSites(seqs->getSequence(name));
    seqNames.push_back(name);
  }


  //Sort matrix entries:
  vector<Index> distances;
  for (unsigned int i = 0; i < dist->size()-1; i++)
    for (unsigned int j = i+1; j < dist->size(); j++)
      distances.push_back(Index((*dist)(i, j), i , j));
  sort(distances.begin(), distances.end());


    ///////////////////////////////////////////////////
    // Performing the actual sequence sampling
    ///////////////////////////////////////////////////
  if (deleteMeth == "random")
  {
    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", phylomerge.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    vector<string> sample(sampleSize);
    RandomTools::getSample(seqNames, sample, false);
    seqNames = sample;
    
    double mini = -log(0.);
    for (unsigned int i =  0; i < seqNames.size() - 1; ++i)
      for (unsigned int j = i + 1; j < seqNames.size(); ++j)
      {
        double d = (*dist)(seqNames[i], seqNames[j]);
        if (d < mini) mini = d;
      }
    ApplicationTools::displayResult("Minimal distance in final data set:", mini);
  }
  else if (deleteMeth == "threshold")
  {
    double threshold = ApplicationTools::getDoubleParameter("threshold", phylomerge.getParams(), 0.01);
    ApplicationTools::displayResult("Distance threshold", threshold);

    unsigned int rm = 0;
    while (distances[0].distance <= threshold)
    {
      //We need to chose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      unsigned int pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + pos); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
      if (distances.size() == 0)
        throw Exception("Error, all sequences have been removed with this criterion!");
    }
    ApplicationTools::displayResult("Number of sequences kept:", seqNames.size());
  }
  else if (deleteMeth == "sample")
  {

    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", phylomerge.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    
    unsigned int rm = 0;
    while (seqNames.size() > sampleSize)
    {
      //We need to choose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      unsigned int pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + pos); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
    }
    ApplicationTools::displayResult("Minimal distance in final data set:", distances[0].distance);
  }
  else if (deleteMeth == "taxon")
  {
    //First we have to build the unrooted tree, in case it was not the input method
    if(inputMethod != "tree") {
      BioNJ *bionj = new BioNJ (*dist, false, true);
		if (tree) delete tree;
      tree = bionj->getTree();
		delete bionj;
    }
	
    //We reroot the tree by a leaf
    Node* newOutgroup = tree->getLeaves()[0];
    tree->newOutGroup(newOutgroup);
    if ( ! tree->isRooted() ) {
        std::cout << "Tree not rooted, rerooting." <<std::endl;
        tree->newOutGroup(tree->getLeaves()[1]);
    }

    std::cout << TreeTemplateTools::treeToParenthesis(*tree) <<std::endl;

	if (! newOutgroup->hasDistanceToFather() ) {
		Node* father = newOutgroup->getFather();
		Node* brother ;
		if (father->getSon(0) == newOutgroup) {
			brother = father->getSon(1);
		}
		else {
			brother = father->getSon(0);
		}
		double dist = brother->getDistanceToFather();
		brother->setDistanceToFather(dist/2);
		newOutgroup->setDistanceToFather(dist/2);
	}	
	else {}


    annotateTreeWithSpecies( *tree, seqSp );
     
    //If we want to rearrange the tree on its unsupported branches
      bool rearrangeTree = ApplicationTools::getBooleanParameter("rearrange.tree", phylomerge.getParams(), false, "", true, false);
 
      if (rearrangeTree) {
          double bootstrapThreshold = ApplicationTools::getDoubleParameter("bootstrap.threshold", phylomerge.getParams(), 0.7);
          if (bootstrapThreshold>1)
          {
              bootstrapThreshold = bootstrapThreshold/100;
          }
          refineTree (*tree,
                      seqSp,
                      speciesToRefine,
                      bootstrapThreshold);
      }
    vector <string> sequencesAncestor;
    seqNames = selectSequencesToKeep(*tree, tree->getRootNode(), sequencesAncestor, critMeth, *seqs, speciesToRefine);

    ApplicationTools::displayResult("Number of sequences kept:", seqNames.size());
  }
  else throw Exception("Unknown deletion method: " + deleteMeth + ".");

  //Write sequences to file:
  AlignedSequenceContainer asc(alphabet);
  for(unsigned int i = 0; i < seqNames.size(); i++)
    if (!asc.hasSequence(seqNames[i])) asc.addSequence(seqs->getSequence(seqNames[i])); //MYSTERY

  SequenceApplicationTools::writeAlignmentFile(asc, phylomerge.getParams());
  
	  if (dist)
		  delete dist;
	  delete tree;
	  delete seqs;
	  delete alphabet;
  phylomerge.done();
  }
  catch (exception& e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}











    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////






/****************************************************************
 * Idea for a new function to cut the tree in subtrees of putative orthologs.
 * Aim: cut at duplication nodes.
 * Algorithm:
 * Double-recursive tree traversal, where each node is annotated
 * with number of genes g in underlying subtree,
 * with number of species s in underlying subtree.
 * Then, additional tree traversal for identifying duplication nodes:
 * N father of n1 and n2 is a duplication node iff:
 * - g(N)-(max(g(n1), g(n2)) > 0 (always true because g(N) = sum (g(n1), g(n2)))
 * - s(N)-max(s(n1), s(n2)) = 0 or (weaker) inter(s(n1), s(n2)) != 0.
 * Then putative duplication nodes have been identified. What to do for each?
 * i) remove subtree that seems least promising (smallest sequences, smallest number of species, largest branches...)
 * ii) assemble all families of putative orthologs, and then think of some reasonable way to select the best family... Can become prohibitively combinatorial.
 * Overall: seems a bit complex to program properly...
 ****************************************************************/


/****************************************************************
 * Idea for a new simple function to cut the tree in subtrees of putative orthologs.
 * - Try to cut at each branch of the tree.
 * - Then compute the size of each subtree, and compute whether it contains at most one sequence per species.
 * -
 ****************************************************************/




  //Outputting a pruned tree, if the sequences have not been merged
    /*TEMP  if ( critMeth != "merge" ) {
          string outTreePath = ApplicationTools::getAFilePath("output.pruned.tree", phylomerge.getParams(), true, false);
        //  DistanceMatrix 
        if (matrix) delete matrix;
        matrix = new DistanceMatrix::DistanceMatrix(seqNames);
          //We build the reduced matrix from the big tree
          const Node * root=tree->getRootNode();
          int *idi, *idj;
          string nomi, nomj;
          for (unsigned int i=0; i<matrix->getNumberOfRows();i++) {
              for (unsigned int j=0; j<matrix->getNumberOfColumns ();j++) {
                  //We compute distances between all the leaves we keep
                  nomi=matrix->getName(i);
                  nomj=matrix->getName(j);   
                  TreeTemplateTools::searchLeaf(*root, nomi, idi);
                  TreeTemplateTools::searchLeaf(*root, nomj, idj);
                  double dij = TreeTools::getDistanceBetweenAnyTwoNodes(*tree, *idi, *idj);
                  matrix->operator()(i,j)=dij;
                  matrix->operator()(j,i)=dij;
              }
          }

          NeighborJoining nj (*matrix, false); //build an unrooted tree
          Newick newick;
          newick.write (*(nj.getTree() ), outTreePath);
      }*/
