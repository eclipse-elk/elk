/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.service;

import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.ui.IWorkbenchPart;

/**
 * Interface for connectors of diagram layout. A diagram layout connector is responsible for transforming
 * the diagram contained in a workbench part into a layout graph, which is an instance of the
 * KGraph meta model. Furthermore it must handle the transfer of concrete layout data from the
 * layout graph back to the diagram after a layout has been computed.
 * 
 * @author msp
 */
public interface IDiagramLayoutConnector {

    /**
     * Build a KGraph instance for the given diagram. The resulting layout graph should reflect
     * the structure of the original diagram. All graph elements must have
     * {@link org.eclipse.elk.core.klayoutdata.KShapeLayout KShapeLayouts} or
     * {@link org.eclipse.elk.core.klayoutdata.KEdgeLayout KEdgeLayouts} attached,
     * and their modification flags must be set to {@code false}.
     * 
     * <p>At least one of the two parameters must be non-null.</p>
     * 
     * <p>This method is usually called from the UI thread.</p>
     * 
     * @param workbenchPart
     *            the workbench part for which layout is performed, or {@code null} if there
     *            is no workbench part for the diagram
     * @param diagramPart
     *            the parent object for which layout is performed, or {@code null} if the
     *            whole diagram shall be layouted
     * @return a layout graph mapping, or {@code null} if the given workbench part or diagram part
     *            is not supported
     */
    LayoutMapping buildLayoutGraph(IWorkbenchPart workbenchPart, Object diagramPart);

    /**
     * Apply the computed layout back to the diagram. Graph elements whose modification flag
     * was not set during layout should be ignored.
     * 
     * <p>This method is usually called from the UI thread.</p>
     * 
     * @param mapping a layout mapping that was created by this layout connector
     * @param settings general settings for applying layout, e.g. whether to use animation
     */
    void applyLayout(LayoutMapping mapping, IPropertyHolder settings);
    
}
