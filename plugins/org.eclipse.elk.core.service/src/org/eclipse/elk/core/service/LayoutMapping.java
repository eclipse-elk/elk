/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.service;

import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.MapPropertyHolder;
import org.eclipse.ui.IWorkbenchPart;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

/**
 * A layout mapping describes the relation between a graphical diagram and the layout graph
 * that is created by diagram layout managers.
 *
 * @author msp
 */
public class LayoutMapping extends MapPropertyHolder {
    
    private static final long serialVersionUID = 2237409212851510612L;
    
    /** the bidirectional mapping of layout graph elements to diagram parts. */
    private final BiMap<ElkGraphElement, Object> graphElemMap = HashBiMap.create();
    /** the top-level parent node of the layout graph. */
    private ElkNode layoutGraph;
    /** the top-level diagram part. */
    private Object parentElement;
    /** the workbench part for wich the mapping was created, if any. */
    private final IWorkbenchPart workbenchPart;
    
    /**
     * Create a layout mapping.
     * 
     * @param theWorkbenchPart the workbench part for which the mapping is created, which may be {@code null}
     */
    public LayoutMapping(final IWorkbenchPart theWorkbenchPart) {
        this.workbenchPart = theWorkbenchPart;
    }
    
    /**
     * Returns the bidirectional mapping of layout graph elements to diagram parts.
     * 
     * @return the graph element map
     */
    public BiMap<ElkGraphElement, Object> getGraphMap() {
        return graphElemMap;
    }
    
    /**
     * Set the top-level parent node of the layout graph.
     * 
     * @param layoutGraph the layout graph
     */
    public void setLayoutGraph(final ElkNode layoutGraph) {
        this.layoutGraph = layoutGraph;
    }
    
    /**
     * Returns the top-level parent node of the layout graph.
     * 
     * @return the layout graph
     */
    public ElkNode getLayoutGraph() {
        return layoutGraph;
    }
    
    /**
     * Set the top-level diagram part.
     * 
     * @param parentElem the parent diagram part
     */
    public void setParentElement(final Object parentElem) {
        this.parentElement = parentElem;
    }
    
    /**
     * Returns the top-level diagram part.
     * 
     * @return the parent diagram part
     */
    public Object getParentElement() {
        return parentElement;
    }
    
    /**
     * Returns the workbench part, or {@code null}.
     */
    public IWorkbenchPart getWorkbenchPart() {
        return workbenchPart;
    }

}
