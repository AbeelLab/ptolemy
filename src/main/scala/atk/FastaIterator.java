package atk;

/**
 * Copied from Thomas Abeel's ATK library on 29/11/2017
 */

import java.io.*;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * This class provides an iterator over a file in FASTA formatting.
 *
 * This format contains a one line header followed by lines of sequence data.
 * Sequences in fasta formatted files are preceded by a line starting with a "
 * &gt;" symbol. The first word on this line is the name of the sequence. The
 * rest of the line is a description of the sequence.
 *
 * After the header line, one or more comments, distinguished by a semi-colon
 * (;) at the beginning of the line, may occur. Most databases and
 * bioinformatics applications do not recognize these comments so their use is
 * discouraged, but they are part of the official format.
 *
 * The remaining lines contain the sequence itself. Blank lines in a FASTA file
 * are ignored, and so are spaces or other gap symbols (dashes, underscores,
 * periods) in a sequence. Fasta files containing multiple sequences are just
 * the same, with one sequence listed right after another. This format is
 * accepted for many multiple sequence alignment programs.
 *
 * @author Thomas Abeel
 *
 */
public class FastaIterator implements Iterable<Record>, Iterator<Record> {

    private BufferedReader in = null;

    private String lastLine = null;

    /**
     * Construct a new fasta iterator based on an inputstream
     *
     * @param is
     *            inputstream
     * @throws IOException
     */
    public FastaIterator(InputStream is) throws IOException {
        in = new BufferedReader(new InputStreamReader(is));
        assert (in != null);
        lastLine = readLine(in);
        next = createNextRecord();
    }

    /**
     * Create a new fasta iterator for the file denoted by the provided name
     *
     * @param file
     *            file name
     * @throws IOException
     */
    public FastaIterator(String file) throws IOException {
        this(new File(file));
    }

    /**
     * Construct a fasta iterator for the file that is provided
     *
     * @param file
     *            file to iterate
     * @throws IOException
     */
    public FastaIterator(File file) throws IOException {
        this(new FileInputStream(file));

    }

    private Record next = null;

    private String readLine(BufferedReader in) throws IOException {
        String line = in.readLine();
        while (line != null && (line.startsWith(";") || line.length() == 0)) {
            line = in.readLine();
        }
        return line;
    }

    private Record createNextRecord() {
        if (lastLine == null)
            return null;
        StringBuffer description = new StringBuffer();
        StringBuffer sequence = new StringBuffer();
        // read description
        if (lastLine.startsWith(">"))
            description.append(lastLine);
        else
            throw new RuntimeException("The fasta file is not valid around the following line:\n\t" + lastLine);
        // read sequence
        try {
            lastLine = readLine(in);
            while (lastLine != null && !lastLine.startsWith(">")) {
                sequence.append(lastLine);
                lastLine = readLine(in);
            }
            if (description != null)
                return new Record(description.toString(), sequence.toString());
            else
                return null;
        } catch (IOException e) {
            return null;
        }
    }

    public Iterator<Record> iterator() {
        return this;
    }

    public boolean hasNext() {
        return next != null;
    }

    public Record next() {
        Record tmp = next;
        try {
            next = createNextRecord();
        } catch (RuntimeException e) {
            throw new NoSuchElementException(e.getMessage());
        }
        return tmp;
    }

    public void remove() {
        throw new UnsupportedOperationException();

    }

    public void close() {
        try {
            in.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

}