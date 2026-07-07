/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.debug.views;

import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.jface.action.IToolBarManager;
import org.eclipse.jface.action.Separator;
import org.eclipse.swt.SWT;
import org.eclipse.swt.custom.ScrolledComposite;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.ui.IViewPart;
import org.eclipse.ui.IWorkbenchPage;
import org.eclipse.ui.IWorkbenchWindow;
import org.eclipse.ui.PlatformUI;
import org.eclipse.ui.part.ViewPart;

/**
 * A viewer for layout graphs.
 * 
 * This class started as a copy of DisCoDebugView.java, commit 71bb8c2f542,
 * 2016-05-07. Changes are commented accordingly.
 */
public class DisCoDebugView extends ViewPart {

    /** the view identifier. */
    public static final String VIEW_ID = "org.eclipse.elk.debug.disCoDebugView"; //$NON-NLS-1$

    /** the scrolled composite that contains the graph canvas. */
    private ScrolledComposite scrolledComposite;
    /**
     * by mic: the canvas used to draw layout graphs and their DCGraph and polyomino
     * representation.
     */
    private DisCoGraphRenderingCanvas graphCanvas;

    /**
     * Creates a layout graph view.
     */
    public DisCoDebugView() {
        super();
    }

    /**
     * Updates the view by painting the given layout graph.
     * 
     * @param graph
     *            the graph to be drawn.
     */
    public static void updateWithGraph(final ElkNode graph) {
        PlatformUI.getWorkbench().getDisplay().syncExec(new Runnable() {
            @Override
            public void run() {
                DisCoDebugView activeView = findView();
                if (graph != null && activeView != null) {
                    ElkNode nodeCopy = EcoreUtil.copy(graph);
                    activeView.getCanvas().setLayoutGraph(nodeCopy);
                }
            }
        });
    }

    /**
     * Tries to find the relevant currently open view.
     */
    private static DisCoDebugView findView() {
        IWorkbenchWindow activeWindow = PlatformUI.getWorkbench().getActiveWorkbenchWindow();
        if (activeWindow == null) {
            return null;
        }

        IWorkbenchPage activePage = activeWindow.getActivePage();
        if (activePage == null) {
            return null;
        }

        // find layout graph view
        IViewPart viewPart = activePage.findView(DisCoDebugView.VIEW_ID);
        if (viewPart instanceof DisCoDebugView) {
            return (DisCoDebugView) viewPart;
        }

        return null;
    }

    @Override
    public void createPartControl(final Composite parent) {
        // create actions in the view toolbar
        IToolBarManager toolBarManager = getViewSite().getActionBars().getToolBarManager();

        // create canvas for layout graphs
        scrolledComposite = new ScrolledComposite(parent, SWT.H_SCROLL | SWT.V_SCROLL);
        // by mic: Another canvas is used instead of GraphRenderingCanvas to not
        // only draw the graph but also its DCGraph and polyomino
        // representation.
        graphCanvas = new DisCoGraphRenderingCanvas(scrolledComposite);
        scrolledComposite.setContent(graphCanvas);
        // by mic: the toolbar allows the user to show/hide the original graph's
        // layout, the associated DCGraph and the lower resolution polyominoes.
        toolBarManager.add(new SelectivePaintToggleAction("elkgraph", "Toggle ElkGraph", "Shows and hides the ElkGraph",
                graphCanvas));
        toolBarManager.add(new SelectivePaintToggleAction("dcgraph", "Toggle DCGraph", "Shows and hides the DCGraph",
                graphCanvas));

        toolBarManager.add(new SelectivePaintToggleAction("polyominoes", "Toggle Polyominoes",
                "Shows and hides polyominoes", graphCanvas));
        toolBarManager.add(new Separator());

        DisCoGraphRenderer renderer = graphCanvas.getRenderer();
        ConfigState state = renderer.getState();
        Integer upperBound = state.getMaxDepth();

        upperBound = upperBound == null ? 0 : upperBound - 1;
        toolBarManager
                .add(new DropdownIntRangeAction(state, "Show levels ", 0, upperBound, 0, State.LOWER, graphCanvas));
        toolBarManager
                .add(new DropdownIntRangeAction(state, "to ", 0, upperBound, upperBound, State.UPPER, graphCanvas));

        toolBarManager.add(new Separator());
        toolBarManager.add(new DisCoImageExportAction(this));

    }

    @Override
    public void setFocus() {
        scrolledComposite.setFocus();
    }

    // by mic: Another canvas is used instead of GraphRenderingCanvas to not
    // only draw the graph but also its DCGraph and polyomino
    // representation.
    /**
     * Returns the DisCo graph canvas.
     * 
     * @return the DisCo graph canvas
     */
    public DisCoGraphRenderingCanvas getCanvas() {
        return graphCanvas;
    }

}
