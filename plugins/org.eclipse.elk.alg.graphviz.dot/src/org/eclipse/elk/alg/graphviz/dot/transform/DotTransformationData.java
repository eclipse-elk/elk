/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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
package org.eclipse.elk.alg.graphviz.dot.transform;

import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.graph.properties.MapPropertyHolder;

/**
 * Pretty much a copy of the format service's {@code TransformationData} class. This class is necessary
 * to keep GraphViz Dot itself independent from the formats service.
 * 
 * @param <S> source graph type.
 * @param <T> target graph type.
 * @author cds
 */
public final class DotTransformationData<S, T> extends MapPropertyHolder
        implements IDotTransformationData<S, T> {

    /** Serialization ID. */
    private static final long serialVersionUID = 2334583523009191492L;
    
    /** the original source graph. */
    private S sourceGraph;
    /** the transformed target graphs. */
    private final List<T> layoutGraphs = new LinkedList<T>();
    /** the log messages. */
    private final List<String> logMessages = new LinkedList<String>();
    

    @Override
    public void setSourceGraph(final S srcGraph) {
        this.sourceGraph = srcGraph;
    }

    @Override
    public S getSourceGraph() {
        return sourceGraph;
    }

    @Override
    public List<T> getTargetGraphs() {
        return layoutGraphs;
    }

    @Override
    public void log(final String msg) {
        logMessages.add(msg);
    }

    @Override
    public Iterable<String> getMessages() {
        return logMessages;
    }
}
