/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.layouter.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * An input stream that forwards all read data to an output stream.
 * 
 * @author msp
 */
public class ForwardingInputStream extends InputStream {

    /** the input stream that is wrapped by this instance. */
    private final InputStream inputStream;
    /** the output stream to which all read data is forwarded. */
    private final OutputStream outputStream;
    
    /**
     * Creates a forwarding input stream with the given input and output streams.
     * 
     * @param theinputStream the input stream from which data is read
     * @param theoutputStream the output stream to which data is written
     */
    public ForwardingInputStream(final InputStream theinputStream,
            final OutputStream theoutputStream) {
        this.inputStream = theinputStream;
        this.outputStream = theoutputStream;
    }
    
    @Override
    public int read() throws IOException {
        int data = inputStream.read();
        if (data >= 0) {
            outputStream.write(data);
        }
        return data;
    }
    
    @Override
    public int read(final byte[] b) throws IOException {
        int count = inputStream.read(b);
        if (count > 0) {
            outputStream.write(b, 0, count);
        }
        return count;
    }
    
    @Override
    public int read(final byte[] b, final int off, final int len) throws IOException {
        int count = inputStream.read(b, off, len);
        if (count > 0) {
            outputStream.write(b, off, count);
        }
        return count;
    }
    
    @Override
    public long skip(final long n) throws IOException {
        return inputStream.skip(n);
    }

    @Override
    public int available() throws IOException {
        return inputStream.available();
    }
    
    /**
     * Closes this input stream and releases any system resources associated
     * with the stream. This method does <em>not</em> close the output stream.
     * 
     * @throws IOException if an I/O error occurs.
     */
    public void close() throws IOException {
        inputStream.close();
    }
    
    @Override
    public synchronized void mark(final int readlimit) {
        inputStream.mark(readlimit);
    }
    
    @Override
    public synchronized void reset() throws IOException {
        inputStream.reset();
    }
    
    @Override
    public boolean markSupported() {
        return inputStream.markSupported();
    }

}
