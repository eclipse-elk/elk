/*******************************************************************************
 * Copyright (c) 2012, 2015 Kiel University and others.
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
package org.eclipse.elk.conn.gmf;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.core.runtime.IAdaptable;
import org.eclipse.core.runtime.Platform;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.gmf.runtime.common.core.service.AbstractProvider;
import org.eclipse.gmf.runtime.common.core.service.IOperation;
import org.eclipse.gmf.runtime.diagram.ui.editparts.AbstractBorderItemEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.IGraphicalEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.ShapeNodeEditPart;
import org.eclipse.gmf.runtime.diagram.ui.render.editparts.AbstractImageEditPart;
import org.eclipse.gmf.runtime.diagram.ui.services.layout.ILayoutNode;
import org.eclipse.gmf.runtime.diagram.ui.services.layout.ILayoutNodeOperation;
import org.eclipse.gmf.runtime.diagram.ui.services.layout.ILayoutNodeProvider;
import org.eclipse.gmf.runtime.notation.Node;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterables;

/**
 * GMF layout provider that executes ELK layout. This enables the execution of ELK layout
 * using the default "Arrange All" button.
 *
 * @author msp
 * @author chsch
 */
public class ElkLayoutProvider extends AbstractProvider implements ILayoutNodeProvider {

    @Override
    public boolean provides(final IOperation operation) {
        return operation instanceof ILayoutNodeOperation;
    }

    /**
     * {@inheritDoc}<br>
     * <br>
     * This method returns <code>true</code> only for {@link IGraphicalEditPart IGraphicalEditParts}
     * whose figures are visible and that contain at least a {@link ShapeNodeEditPart}, which is not
     * an {@link AbstractBorderItemEditPart} or an {@link AbstractImageEditPart}.<br>
     * <br>
     * In short, returns only true if the provided edit part contains children that can be layouted.
     */
    @SuppressWarnings("rawtypes")
    public boolean canLayoutNodes(final List layoutNodes, final boolean shouldOffsetFromBoundingBox,
            final IAdaptable layoutHint) {
        Object o = layoutHint.getAdapter(IGraphicalEditPart.class);

        if (!(o instanceof IGraphicalEditPart)) {
            return false;
        }
        final IGraphicalEditPart parent = (IGraphicalEditPart) o;

        // check the availability of children that can be arranged
        //  computing layout on the element wouldn't make sense otherwise
        if (layoutNodes.isEmpty()) {
            return Iterables.any((List<?>) parent.getChildren(), new Predicate<Object>() {
                public boolean apply(final Object o) {
                    return o instanceof ShapeNodeEditPart
                            && !(o instanceof AbstractBorderItemEditPart)
                            && !(o instanceof AbstractImageEditPart);
                }
            });
        } else {
            return Iterables.any((List<?>) layoutNodes, new Predicate<Object>() {
                public boolean apply(final Object o) {
                    ILayoutNode layoutNode = (ILayoutNode) o;
                    IGraphicalEditPart editPart = findEditPart(parent, layoutNode.getNode());
                    return editPart instanceof ShapeNodeEditPart
                            && !(editPart instanceof AbstractBorderItemEditPart)
                            && !(editPart instanceof AbstractImageEditPart);
                }
            });
        }
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Runnable layoutLayoutNodes(final List layoutNodes, final boolean offsetFromBoundingBox,
            final IAdaptable layoutHint) {
        // fetch general settings from preferences
        final boolean zoomToFit = Platform.getPreferencesService().getBoolean(
                "org.eclipse.elk.core.ui",
                "org.eclipse.elk.zoomToFit", false, null);
        final boolean progressDialog = Platform.getPreferencesService().getBoolean(
                "org.eclipse.elk.core.ui",
                "org.eclipse.elk.progressDialog", false, null);
        
        // determine the elements to process
        final Object diagramPart;
        if (layoutNodes.isEmpty()) {
            diagramPart = layoutHint.getAdapter(IGraphicalEditPart.class);
        } else {
            List<IGraphicalEditPart> partList = new ArrayList<IGraphicalEditPart>(layoutNodes.size());
            IGraphicalEditPart parent = (IGraphicalEditPart) layoutHint.getAdapter(
                    IGraphicalEditPart.class);
            if (parent != null) {
                Iterator<?> nodeIter = layoutNodes.iterator();
                while (nodeIter.hasNext()) {
                    ILayoutNode layoutNode = (ILayoutNode) nodeIter.next();
                    IGraphicalEditPart editPart = findEditPart(parent, layoutNode.getNode());
                    if (editPart != null) {
                        partList.add(editPart);
                    }
                }
            }
            diagramPart = partList;
        }
        
        return new Runnable() {
            public void run() {
                DiagramLayoutEngine.invokeLayout(null, diagramPart, false, progressDialog,
                        false, zoomToFit);
            }
        };
    }
    
    /**
     * Find an edit part with the given notation node.
     * 
     * @param parent the parent edit part to start the search
     * @param notationNode a notation node
     * @return the corresponding edit part
     */
    @SuppressWarnings("unchecked")
    private static IGraphicalEditPart findEditPart(final IGraphicalEditPart parent,
            final Node notationNode) {
        LinkedList<IGraphicalEditPart> editPartQueue = new LinkedList<IGraphicalEditPart>();
        editPartQueue.add(parent);
        do {
            IGraphicalEditPart editPart = editPartQueue.removeFirst();
            if (notationNode.equals(editPart.getNotationView())) {
                return editPart;
            }
            editPartQueue.addAll((Collection<? extends IGraphicalEditPart>) editPart.getChildren());
        } while (!editPartQueue.isEmpty());
        return null;
    }
    
}
