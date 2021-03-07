# include <iostream>
# include <fstream>
# include <map>
# include <vector>
# include <string>
# include <queue>
# include <functional>
# include <algorithm>
# include "Wrapper.h"

class ibstream{
    public:
        ibstream( istream & is );
        int readBit ( ) ;
        istream & getInputStream( ) const;
    private:
        istream & in; // The underlying input stream
        char buffer; // Buffer to store eight bits at a time
        int bufferPos; // Position in buffer for next read
};

class obstream
 {
    public:
        obstream( ostream & os ) ;
        ~obstream ( ) ;

        void writeBit( int val ) ;
        void writeBits( const vector<int> & val ) ;
        void flush ( ) ;
        ostream & getOutputStream( ) const;

    private:
        ostream & out; // The underlying output stream
        char buffer; // Buffer to store eight bits at a time
        int bufferPos; // Position in buffer for next write
 };

 static const int BITS_PER_CHAR = 8;
 static const int DIFF_CHARS = 256;

// return bit at position pos in a packed set of bits (pack).

int getBit(char pack, int pos){
    return (pack & (1<<pos)) ? 1 : 0;
}

// Set bit at position pos in a packed set of bits (pack).
void setBit(char & pack, int pos, int val){
    if (val ==1)
        pack |=(val <<pos);
}

// Construct the bit-input stream.
ibstream::ibstream (istream &is):bufferPos(BITS_PER_CHAR), in(is){

}
//read one bit
int ibstream::readBit(){
    if (bufferPos == BITS_PER_CHAR){
        in.get( buffer );
        if (in.eof())
            return EOF;
        bufferPos = 0;
    }
    return getBit(buffer, bufferPos++);
}

// Return underlying input stream.

istream & ibstream::getInputStream() const{
    return in;
}
// Construct the bit-output stream.
obstream::obstream (ostream &os):bufferPos(0), buffer(0), out(os){

}

// Destructor (writes any bits left in the buffer).

obstream::~obstream(){
    flush();
}


// If the buffer is non-empty, write one more char

void obstream::flush(){
    if (bufferPos == 0)
        return;
    out.put(buffer);
    bufferPos = 0;
    buffer = 0;

}
//write one bit
void obstream::writeBit(int val){
    setBit(buffer, bufferPos++, val);
    if (bufferPos == BITS_PER_CHAR)
        flush();
}

//write a vector of bits 
void obstream::writeBits(const vector<int> & val){
    for (int i = 0; i < val.size(); i++)
        writeBit(val[i]);
}

//return underlying output stream

ostream & obstream::getOutputStream() const{
    return out;
}


class CharCounter{
    public:
        CharCounter();
        CharCounter(istream & input);

        int getCount(char ch) const;
        void setCount(char ch, int count);

    private:
        map<char, int, less<char> > theCounts;
};

// Constructor: All counts are zero.

CharCounter::CharCounter()/*: in(cin)*/{}

// Constructor: Get counts by reading from input stream.
CharCounter::CharCounter(istream & input){
    char ch;
    while(!input.get(ch).eof())
        theCounts[ch]++;
}

//return the character count for ch

int CharCounter::getCount(char ch) const{
    map<char,int,less<char>>::const_iterator itr;
    itr = theCounts.find(ch);
    return itr !=theCounts.end() ? (*itr).second : 0;
}

//set the character count for char

void CharCounter::setCount(char ch, int count){
    theCounts[ch] = count;
}

struct HuffNode
{
    int value;
    int weight;
    HuffNode* left;
    HuffNode* right;
    HuffNode* parent;

    HuffNode(int v, int w, HuffNode *lt, HuffNode *rt, HuffNode *pt): value(v), weight(w), left(lt), right(rt), parent(pt){}
};


class HuffmanTree{
    public:
        HuffmanTree();
        HuffmanTree(const CharCounter & cc);

        enum {ERROR = -3, INCOMPLETE = -2, END = -1};

        // Here, vector<int> is usable by ibstream and obstreams

        vector<int> getCode(int ch) const;
        int getChar (const vector<int> & code ) const;

        // Write the encoding table using character counts.

        void writeEncodingTable(ostream & out);
        void readEncodingTable(istream & in);

    private:
        CharCounter theCounts;
        map<int, HuffNode *, less<int>> theNodes;
        HuffNode* root;

        void createTree();
};



// Construct the tree given a CharCounter objecr.
 // Tree will be usable.
HuffmanTree::HuffmanTree( const CharCounter & cc ): theCounts( cc )
 {
 root = NULL;
 createTree ( ) ;
 }

 // Construct the tree in an unusable state.
 // A call to readEncodingTable is expected to follow.
 HuffmanTree::HuffmanTree( )
 {
     root = NULL;
 }
 


 // Return the code corresponding to character ch.
 // (The parameter is an int to accommodate EOF).
 // If code is not found, return a vector of size 0.
 vector<int> HuffmanTree::getCode( int ch ) const
 {
    map<int,HuffNode *,less<int> >::const_iterator itr;
    itr = theNodes.find( ch );
    if( itr == theNodes.end( ) )
        return vector<int>( );
    HuffNode *current = (*itr) .second;

    vector<int> v;
    HuffNode *par = current->parent;
    while( par != NULL )
    {
        if( par->left == current )
            v.push_back( 0 ); // current is a left child
        else
            v.push_back( 1 ); // current is a right child
        current = current->parent;
        par = current->parent;
    }

 reverse( v.begin( ), v.end( ) );
 return v;
 }

 int HuffmanTree::getChar(const vector<int> & code) const{
     HuffNode *p = root;
     for (int i =0; p!= NULL && i < code.size(); i++){
         if (code[i] == 0)
            p=p->left;
        else
            p = p->right;
     }
     if (p==NULL)
        return ERROR;

    return p->value;
 }

 // Write an encoding table to an output stream.
 // Format is character, count (formatted), newline.
 // A zero count terminates the encoding table.

 void HuffmanTree::writeEncodingTable(ostream & out){
     for (int i =0; i<DIFF_CHARS; i++)
        if (theCounts.getCount(i)>0)
            out << static_cast<char>(i)<<theCounts.getCount(i)<<'\n';
    out <<'\0' << 0<<'\n';
 }

 // Read the encoding table from an input stream in format
// given above and then construct the Huffman tree.
// Stream will then be positioned to read compressed data.

void HuffmanTree::readEncodingTable(istream & in){
    for (int i =0; i < DIFF_CHARS; i++)
        theCounts.setCount(i,0);

    char ch, nl;

    int num;

    for (;;){
        in.get(ch);
        in>>num;
        in.get(nl);
        if (num == 0)
            break;
        theCounts.setCount(ch,num);
    }
    createTree();
}

// Comparison function for HuffNode.
// Meaning is reversed so priority-queue will retrive the min.

bool operator<(const HuffNode & lhs, const HuffNode & rhs){
    return lhs.weight > rhs.weight;
}

//construct the huffman coding tree

void HuffmanTree::createTree(){
    priority_queue<Pointer<HuffNode>, vector<Pointer<HuffNode> >, less<Pointer<HuffNode>>> pq;
    for (int i =0; i < DIFF_CHARS; i++)
        if (theCounts.getCount(i) > 0){
            HuffNode *newNode = new HuffNode(i,theCounts.getCount(i), NULL, NULL, NULL);
            theNodes[i] = newNode;
            pq.push(Pointer<HuffNode>(newNode));
        }

    theNodes[END] = new HuffNode(END,1,NULL,NULL,NULL);
    pq.push(Pointer<HuffNode>(theNodes[END]));

    while(pq.size()>1){
        HuffNode *n1 = pq.top(); pq.pop();
        HuffNode *n2 = pq.top(); pq.pop();
        HuffNode *result = new HuffNode(INCOMPLETE, n1->weight + n2->weight,n1,n2,NULL);
        n1->parent = n2->parent = result;
        pq.push(Pointer<HuffNode>(result));

    }
    root = pq.top();
}

class Compressor {
    public:
        static void compress(const string & inFile);
        static void uncompress(const string & compressedFile);

    private:
        static const /*int*/ std::ios_base::openmode READ_MODE;
        static const /*int*/ std::ios_base::openmode WRITE_MODE;
};


const /*int*/ std::ios_base::openmode Compressor::READ_MODE = ios::in | ios::binary;
const /*int*/ std::ios_base::openmode Compressor::WRITE_MODE = ios::out | ios::binary;

// Compress inFile; writes result to a file whose name is
// formed by appending ".hufN.
// Very little error checking is performed.

void Compressor::compress(const string & inFile){
    string compressedFile = inFile + ".huf";
    ifstream in(inFile.c_str(), READ_MODE);

    CharCounter countObj(in);
    HuffmanTree codeTree (countObj);

    ofstream out(compressedFile.c_str(),WRITE_MODE);
    codeTree.writeEncodingTable(out);

    in.clear( ); in.seekg( 0, ios::beg );
    obstream bout ( out );

    char ch;
    while(in.get(ch))
        bout.writeBits(codeTree.getCode(ch & 0xff));
    bout.writeBits(codeTree.getCode(EOF));
}

// Uncompress a file. Write the result to a file whose name
// is formed by adding ".ucU (in reality we would simply
// form the new name by stripping off ".hufN.
void Compressor::uncompress(const string & compressedFile){
    int i;
    string inFile, extension;

    for (i = 0; i < compressedFile.length() - 4; i++)
        inFile += compressedFile[i];
    for (; i<compressedFile.length(); i++)
        extension += compressedFile[i];
    if (extension != ".huf")
    {
        cerr<<"Not a compressed file"<<endl;
        return;
    }

    inFile +=".uc";

    ifstream in(compressedFile.c_str(), READ_MODE);
    ofstream out(inFile.c_str(), WRITE_MODE);

    HuffmanTree codeTree;
    codeTree.readEncodingTable(in);

    ibstream bin(in);
    vector<int> bits;
    int bit;
    int decode;

    for (;;){
        bit = bin.readBit();
        bits.push_back(bit);

        decode = codeTree.getChar(bits);
        if(decode == HuffmanTree::INCOMPLETE)
            continue;
        else if(decode == HuffmanTree::ERROR)
        {
            cerr<<"Error Decoding"<<endl;
            break;
        }
        else if (decode == HuffmanTree::END)
            break;
        else{
            out.put(static_cast<char>(decode));
            bits.resize(0);
        }
    }
}

int main(int argc, char*argv[]){
    if (argc<3){
        cerr<<"usage:" << argv[0]<<"-[cu] files" << endl;
        return 1;
    }
    string option = argv[1];
    for (int i =2; i < argc; i++){
        string nextFile = argv[i];
        if (option == "-c")
            Compressor::compress(nextFile);
        else if (option == "-u")
            Compressor::uncompress(nextFile);

        else
        {
            cerr<<"Illegal option; usage:"<<argv[0]<<" -[cu] files" <<endl;
            return 1;
        }
    }
    return 0;
}

