/*******************************************************************************
 * Copyright (c) 2013, 2022 Kiel University and others.
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
package org.eclipse.elk.alg.libavoid.server;

import java.util.LinkedList;

import org.eclipse.elk.alg.libavoid.server.LibavoidServer.Cleanup;


/**
 * A pool for Libavoid server process instances.
 *
 * @author msp
 */
public final class LibavoidServerPool {
    
    /** the singleton instance of the server pool. */
    public static final LibavoidServerPool INSTANCE = new LibavoidServerPool();
    
    /**
     * Hide constructor to avoid instantiation from outside.
     */
    private LibavoidServerPool() {
    }
    
    /** the list of currently available servers. */
    private LinkedList<LibavoidServer> servers = new LinkedList<LibavoidServer>();
    
    /**
     * Fetch an Libavoid server process from the pool, creating one if necessary.
     * 
     * @return an Libavoid server process
     */
    public LibavoidServer fetch() {
        synchronized (servers) {
            if (servers.isEmpty()) {
                return new LibavoidServer();
            }
            return servers.removeFirst();
        }
    }
    
    /**
     * Release a previously created server process into the pool. Only instances that
     * are still usable may be released. Process instances that lead to errors must be closed
     * without releasing them.
     * 
     * @param server an Libavoid server process
     */
    public void release(final LibavoidServer server) {
        synchronized (servers) {
            servers.add(server);
        }
    }
    
    /**
     * Dispose all created server instances.
     */
    public void dispose() {
        synchronized (servers) {
            for (LibavoidServer server : servers) {
                server.cleanup(Cleanup.STOP);
            }
            servers.clear();
        }
    }

}
