//
// File: PhyloCut.cpp
// Created by: Julien Dutheil, Bastien Boussau
// Created on: Sunday, December 2nd 2007 16:48
//

/*
Copyright or  or Copr. CNRS

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
//#include <../../home/boussau/Programs/Boost/boost_1_57_0/boost/concept_check.hpp>

const
std::string
ONE = "E";						 //Property used to store the taxon corresponding to parent and son0 nodes
const
std::string
TWO = "Fu";						 //Property used to store the taxon corresponding to parent and son1 nodes
const
std::string
THREE = "S";					 //Property used to store the taxon corresponding to son0 and son1 nodes

//Compilation:
//g++  -pipe -o phylocut PhyloCut.cpp -I/usr/local/include  -L. -L/usr/local/lib  -g -Wall -fopenmp -std=c++0x -lbpp-core -lbpp-seq -lbpp-phyl

using
namespace
bpp;

void
help ()
{
	(*ApplicationTools::message <<
		"__________________________________________________________________________").endLine
		();
	(*ApplicationTools::message <<
		"phylocut parameter1_name=parameter1_value").endLine ();
	(*ApplicationTools::message <<
		"      parameter2_name=parameter2_value ... param=option_file").endLine ();
	(*ApplicationTools::message).endLine ();
	(*ApplicationTools::message << "      Options considered: ").endLine ();
	(*ApplicationTools::message <<
		" - input.tree.file=file. File with the tree in it, format NHX.").endLine ();
  (*ApplicationTools::message <<
    " - input.sequence.file=file. File with the sequence alignment in it. Output sequence files will be written to files with the same base name.").endLine ();
	(*ApplicationTools::message <<
		" - ids.to.cut=file. File giving the list of species ids for cutting. If not provided, all ids are considered. If an id is included in the list of ids, then all duplication nodes with this species id will be used for cutting, meaning that the two subtrees of duplication nodes of the said id will be put in independent alignments, and independent trees.").endLine
		();
	(*ApplicationTools::message <<
		"__________________________________________________________________________").endLine
		();
}


/**************************************************************************/
/*Removes a leaf in a tree*/
void removeLeaf ( TreeTemplate<Node> & tree, std::string toRemove ) {
  // std::cout <<"in removeLeaf 1"<<std::endl;
  Node * NToRemove  = tree.getNode ( toRemove );
  if ( !NToRemove->hasFather() ) {
    // std::cout <<"node is root !!!!!"<<std::endl;
    Node * father = NToRemove->getSon ( 0 );
    Node * son1 = father->getSon ( 1 );
    father->removeFather();
    tree.newOutGroup ( son1->getId() );
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
  // std::cout <<"in removeLeaf 2"<<std::endl;
  Node * father = NToRemove->getFather();
  // std::cout <<"in removeLeaf 3"<<std::endl;
  Node * brother;
  for ( unsigned int i=0; i<father->getNumberOfSons(); i++ )
    if ( father->getSon ( i ) !=NToRemove ) {
      brother=father->getSon ( i );
      break;
    }
    // std::cout <<"in removeLeaf 4"<<std::endl;
    double distBro;
  try {
    distBro = brother->getDistanceToFather();
  }
  catch ( std::exception ) {
    distBro = 0.000001;
  }
  // std::cout <<"in removeLeaf 5"<<std::endl;
  if ( !father->hasFather() ) {
    // std::cout <<"father is root !!!!!"<<std::endl;
    // brother->removeFather();
    tree.rootAt ( brother->getId() );
    // std::cout <<"After rootAt"<<std::endl;
    //tree.newOutGroup(brother->getSon(0)->getId());
    for ( unsigned int i = 0; i<brother->getNumberOfSons(); i++ ) {
      if ( brother->getSon ( i ) ==NToRemove ) {
        brother->removeSon ( i );
        break;
      }
    }
    tree.newOutGroup ( brother->getSon ( 0 )->getId() );
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
  double distFa;
  try {
    distFa = father->getDistanceToFather();
  }
  catch ( std::exception ) {
    distFa = 0.000001;
  }
  // std::cout <<"in removeLeaf 6"<<std::endl;
  Node * grandFather = father->getFather();
  grandFather->addSon ( brother );
  brother->setDistanceToFather ( distBro+distFa );
  for ( unsigned int i=0; i<grandFather->getNumberOfSons(); i++ )
    if ( grandFather->getSon ( i ) ==father ) {
      grandFather->removeSon ( i );
      break;
    }
    /*if (!grandFather->hasFather()) {
     *      tree.newOutGroup(grandFather->getSon(0)->getId());
}*/
    tree.resetNodesId();
    /*  delete NToRemove;
     *      delete father;*/
}

/**************************************************************************/


TreeTemplate<Node>* extractSubtree( TreeTemplate<Node>* tree, std::vector<std::string> names) {
  TreeTemplate<Node>* newTree = tree->clone();
  std::vector<std::string> allNames = newTree->getLeavesNames();
  std::vector<std::string> namesToRemove ;
  VectorTools::diff(allNames, names, namesToRemove);
 for (size_t i = 0; i<namesToRemove.size(); ++i) {
 //  std::cout << " names i "<< namesToRemove[i] << " and size: "<< newTree->getNumberOfLeaves() <<std::endl;
   if (newTree->getNumberOfLeaves() > 2)
      removeLeaf ( *newTree, namesToRemove[i] );
   else 
     TreeTemplateTools::dropLeaf( *newTree, namesToRemove[i] );

//      TreeTemplateTools::dropLeaf( *newTree, namesToRemove[i] );

 }
  return newTree;
}


void  cutAtDuplicationNodes(Node* node, 
                            VectorSiteContainer *seqs, 
                            bool cutEverywhere, 
                            std::set<int> spIdsToCut, 
                            std::vector<VectorSiteContainer*>& subAlns
                           ) {

      if ( (node->hasBranchProperty("Ev") && (dynamic_cast<const BppString*>(node->getBranchProperty("Ev") ) )->toSTL() == "D") || ( (node->hasBranchProperty("D") && (dynamic_cast<const BppString*>(node->getBranchProperty("D") ) )->toSTL() == "Y") ) ) {
          bool weCut = false;
          if ( cutEverywhere ) {
            weCut = true;
          }
        else {
          int spId =  TextTools::toInt ( (dynamic_cast<const BppString*>(node->getNodeProperty( "S" ) ) )->toSTL() );
          if ( spIdsToCut.find(spId) != spIdsToCut.end() ) {
            weCut = true;
          }
        }

        if (weCut) {

          //Create an alignment for each subtree
          Node* son0 = node->getSon (0);
          Node* son1 = node->getSon (1);
          size_t numAln0 = subAlns.size();
          size_t numAln1 = subAlns.size()+1;
          std::vector<std::string> leafNames0 = TreeTemplateTools::getLeavesNames( *son0 );
          VectorSiteContainer* aln0 = new VectorSiteContainer(seqs->getAlphabet());
          std::vector<std::string> leafNames1 = TreeTemplateTools::getLeavesNames( *son1 );

          for (size_t i = 0; i< leafNames0.size() ; ++i) {
              aln0->addSequence( *(seqs->removeSequence(leafNames0[i] ) ) );
          }
          VectorSiteContainer* aln1 = new VectorSiteContainer(seqs->getAlphabet());
          for (size_t i = 0; i< leafNames1.size() ; ++i) {
              aln1->addSequence( *(seqs->removeSequence(leafNames1[i]) ) );
          }
          subAlns.push_back( aln0 );
          subAlns.push_back( aln1 );
          cutAtDuplicationNodes (son0, subAlns[numAln0], cutEverywhere, spIdsToCut, subAlns);
          cutAtDuplicationNodes (son1, subAlns[numAln1], cutEverywhere, spIdsToCut, subAlns);

       }
       else {
         if ( ! node->isLeaf() ) {
            Node* son0 = node->getSon (0);
            Node* son1 = node->getSon (1);
            cutAtDuplicationNodes (son0, seqs, cutEverywhere, spIdsToCut, subAlns);
            cutAtDuplicationNodes (son1, seqs, cutEverywhere, spIdsToCut, subAlns);
        }
       }
      }
      else {

        if ( ! node->isLeaf() ) {
            Node* son0 = node->getSon (0);
            Node* son1 = node->getSon (1);
            cutAtDuplicationNodes (son0, seqs, cutEverywhere, spIdsToCut, subAlns);
            cutAtDuplicationNodes (son1, seqs, cutEverywhere, spIdsToCut, subAlns);
        }
      }
    return;
   
}




void  cutAtDuplicationNodes(TreeTemplate<Node>* tree, 
                            VectorSiteContainer *seqs, 
                            bool cutEverywhere, 
                            std::set<int> spIdsToCut, 
                            std::string seqName, 
                            std::string  treeName
                           ) {

   Node* root = tree->getRootNode();
   std::vector<VectorSiteContainer*> subAlns;
   subAlns.push_back(seqs);

   cutAtDuplicationNodes(root, subAlns[0], cutEverywhere, spIdsToCut, subAlns);
   //Now we have a list of subalignments that we need to write to file, 
   //along with their associated trees.
   string sequenceFormat =  "Fasta";
   BppOAlignmentWriterFormat bppoWriter( 0 );
   auto_ptr<OAlignment> oAln(bppoWriter.read(sequenceFormat));

   Nhx* nhx = new Nhx();
   std::ofstream out;
   if (subAlns.size() > 0) {
    for (int i = subAlns.size()-1; i>=0; i--) {
        if (subAlns[i]->getNumberOfSequences() != 0 ) {
          std::string sequenceFilePath = seqName + TextTools::toString(i);

          // Write sequences:
          oAln->writeAlignment(sequenceFilePath, *(subAlns[i]), true);
          std::vector<string> names = subAlns[i]->getSequencesNames();
          TreeTemplate<Node>* newTree = nullptr;
          if (tree->getNumberOfLeaves() > names.size() ){
            newTree = extractSubtree(tree, names);
          }
          else {

            newTree = tree->clone();
          }
          
          if (newTree->getNumberOfLeaves() > 1 && newTree->getNumberOfLeaves() < tree->getNumberOfLeaves()) {
          std::string treeFilePath = treeName + TextTools::toString(i);
          out.open (treeFilePath.c_str(), std::ios::out);
          nhx->write(*newTree, out);
          out.close();
          }
          else {
                   //  std::cout << "Do nothing." <<std::endl; 
          }
        }
        else{
        // std::cout << "No sequence in the alignment." <<std::endl; 
        }

    }
  }
   
}




int
main (int args, char **argv)
{
	cout << "******************************************************************"
		<< endl;
	cout << "*        Bio++ Phylogenetic Cutter, version 0.1                  *"
		<< endl;
	cout << "* Author: B. Boussau                        Last Modif. 08/12/14 *"
		<< endl;
	cout << "******************************************************************"
		<< endl;
	cout << endl;

	if (args == 1)
	{
		help ();
		return 0;
	}

	try
	{

		BppApplication phylocut (args, argv, "PhyloCut");
		phylocut.startTimer ();

		//Get sequences:
		Alphabet *
			alphabet =
			SequenceApplicationTools::getAlphabet (phylocut.getParams ());
		VectorSiteContainer *
			seqs = SequenceApplicationTools::getSiteContainer (alphabet,
			phylocut.
			getParams ());
		ApplicationTools::displayResult ("Initial number of sequences:",
			seqs->getNumberOfSequences ());


		TreeTemplate < Node > *tree = 0;
			string
				treePath =
				ApplicationTools::getAFilePath ("input.tree.file",
				phylocut.getParams (), true,
				true);

			ifstream file_stream (treePath.c_str ());
			vector < string > trees;
			int
				tree_i = 0;
								 //  ########## read trees ############
			if (file_stream.is_open ())
			{
				while (!file_stream.eof ())
				{
					string line;
					getline (file_stream, line);
					if (line.find ("(") != line.npos)
					{
						tree_i++;
						trees.push_back (line);
					}
				}
			}
			if (tree_i > 1)
			{
				std::cout <<
					"More than 1 tree in the input tree file, only considering the first tree."
					<< std::endl;
			}

			Nhx* nhx = new Nhx();
			tree = nhx->parenthesisToTree (trees[0]);

		string idFile;
    std::set<int> spIdsToCut;
    bool cutEverywhere = false;
    std::string speciesFile;
			if (ApplicationTools::getAFilePath
				("ids.to.cut", phylocut.getParams (), false,
				true) != "none")
			{
				speciesFile =
					ApplicationTools::getAFilePath ("ids.to.cut",
					phylocut.getParams (), true,
					true);
				//In this file, the format is expected to be as follows :
				/*
				   id1
				   id2
				   id3
				   ...
				 */
				std::ifstream inSpSeq (speciesFile.c_str ());
				std::string line;
        
				while (getline (inSpSeq, line))
				{
          string temp = TextTools::removeSurroundingWhiteSpaces(line);
          if (TextTools::isDecimalInteger(temp) ) {  
            spIdsToCut.insert( TextTools::toInt( temp ) );
          }
          else {
              std::cout << "WARNING: Not recognizing the line "+ temp + " in file "+speciesFile + " as an integer. Ignoring it." <<std::endl;
          }
				}
			}
			else
			{
        cutEverywhere = true;
        std::cout << "Cutting at all duplication nodes." <<std::endl;
			}

			//get the seq name
			std::string seqName = ApplicationTools::getAFilePath ("input.sequence.file",
        phylocut.getParams (), true,
        true);
      //get the tree file name
      std::string treeName = ApplicationTools::getAFilePath ("input.tree.file",
        phylocut.getParams (), true,
        true);

			//Now we do the actual cutting.
			cutAtDuplicationNodes(tree, seqs, cutEverywhere, spIdsToCut, seqName, treeName);

			vector < string > seqNames = seqs->getSequencesNames ();	
		

		delete tree;
		delete seqs;
		delete alphabet;
		phylocut.done ();
	}
	catch (exception & e)
	{
		cout << endl;
		cout << "_____________________________________________________" << endl;
		cout << "ERROR!!!" << endl;
		cout << e.what () << endl;
		return 1;
	}

	return 0;
}

