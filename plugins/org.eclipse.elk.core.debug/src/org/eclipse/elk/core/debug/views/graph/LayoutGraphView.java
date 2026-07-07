/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.graph;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.elk.core.debug.views.AbstractLayoutDebugView;
import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.elk.core.ui.rendering.GraphRenderingCanvas;
import org.eclipse.elk.core.util.LoggedGraph;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.jface.action.IMenuManager;
import org.eclipse.jface.action.IToolBarManager;
import org.eclipse.jface.action.Separator;
import org.eclipse.jface.resource.JFaceResources;
import org.eclipse.jface.viewers.ArrayContentProvider;
import org.eclipse.jface.viewers.DelegatingStyledCellLabelProvider;
import org.eclipse.jface.viewers.ISelectionChangedListener;
import org.eclipse.jface.viewers.IStructuredSelection;
import org.eclipse.jface.viewers.SelectionChangedEvent;
import org.eclipse.jface.viewers.StructuredSelection;
import org.eclipse.jface.viewers.TableViewer;
import org.eclipse.jface.viewers.TreeViewer;
import org.eclipse.swt.SWT;
import org.eclipse.swt.custom.ScrolledComposite;
import org.eclipse.swt.custom.StackLayout;
import org.eclipse.swt.dnd.DND;
import org.eclipse.swt.dnd.FileTransfer;
import org.eclipse.swt.dnd.Transfer;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Text;

/**
 * A viewer for (intermediate) graphs.
 */
public class LayoutGraphView extends AbstractLayoutDebugView {

    /** The view identifier. */
    public static final String VIEW_ID = "org.eclipse.elk.debug.graphView"; //$NON-NLS-1$
    
    // Actions
    private final LayoutUponLoadSettingAction layoutUponLoadSettingAction = new LayoutUponLoadSettingAction();
    private final LoadGraphAction loadGraphAction = new LoadGraphAction(this);
    private final ReloadFromFileAction reloadFromFileAction = new ReloadFromFileAction(this);
    private final ImageExportAction saveImageAction = new ImageExportAction(this);
    
    // UI Controls
    private GraphListLabelProvider graphListLabelProvider = new GraphListLabelProvider();
    private TableViewer graphListViewer;
    private Composite canvasStackComposite;
    private StackLayout canvasStackLayout;
    private ScrolledComposite visualGraphCanvasComposite;
    private GraphRenderingCanvas visualGraphCanvas;
    private Text textualGraphCanvas;
    
    /**
     * Creates a layout graph view.
     */
    public LayoutGraphView() {
        super(VIEW_ID, DiagramLayoutEngine.PREF_DEBUG_LOGGING);
    }
    
    @Override
    public void dispose() {
        super.dispose();
        
        graphListLabelProvider.dispose();
    }
    
    /**
     * Returns the graph rendering canvas if it currently shows anything meaningful.
     */
    public GraphRenderingCanvas getCanvas() {
        if (visualGraphCanvas.getLayoutGraph() != null && canvasStackLayout.topControl == visualGraphCanvasComposite) {
            return visualGraphCanvas;
        } else {
            return null;
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Execution Tree Handling
    
    @Override
    protected void treeSelectionChanged() {
        super.treeSelectionChanged();
        
        // Retrieve all logged graphs from all selected elements
        List<LoggedGraph> loggedGraphs = getSelectedExecutionInfos().stream()
            .flatMap(info -> info.getLoggedGraphs().stream())
            .collect(Collectors.toList());
        
        graphListViewer.setInput(loggedGraphs);
        
        if (!loggedGraphs.isEmpty()) {
            graphListViewer.setSelection(new StructuredSelection(loggedGraphs.get(0)), true);
        }
        
        updateActionEnablement();
    }

    @Override
    protected Predicate<ExecutionInfo> getTreeFilter() {
        return info -> info.hasLoggedGraphs() || info.hasDescendantsWithLoggedGraphs();
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Graph List Handling
    
    /**
     * Called when the selection is changed on our list of logged graphs.
     */
    private void listSelectionChanged() {
        IStructuredSelection selection = graphListViewer.getStructuredSelection();
        if (selection.isEmpty()) {
            visualGraphCanvas.setLayoutGraph(null);
            textualGraphCanvas.setText("");
        } else {
            // We only put LoggedGraph instances into the list, so this cast should never fail
            showLoggedGraph((LoggedGraph) selection.getFirstElement());
        }
        
        updateActionEnablement();
    }

    /**
     * Fills the graph view depending on the graph's type.
     */
    private void showLoggedGraph(LoggedGraph graphInfo) {
        switch (graphInfo.getGraphType()) {        
        case ELK:
            showGraphOnCanvas((ElkNode) graphInfo.getGraph());
            break;
            
        default:
            showTextOnCanvas(graphInfo.getGraph().toString());
            break;
        }
    }

    private void showGraphOnCanvas(ElkNode graph) {
        // putting visual graph view in front and set graph
        canvasStackLayout.topControl = visualGraphCanvasComposite;
        canvasStackComposite.layout();
        visualGraphCanvas.setLayoutGraph(graph);
        textualGraphCanvas.setText("");
    }

    private void showTextOnCanvas(String text) {
        canvasStackLayout.topControl = textualGraphCanvas;
        canvasStackComposite.layout();
        textualGraphCanvas.setText(text);
        visualGraphCanvas.setLayoutGraph(null);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // UI Setup
    
    @Override
    protected void customizeMenu(IMenuManager menuManager) {
        menuManager.add(layoutUponLoadSettingAction);
    }
    
    @Override
    protected void customizeToolBar(IToolBarManager toolBarManager) {
        toolBarManager.add(loadGraphAction);
        toolBarManager.add(reloadFromFileAction);
        toolBarManager.add(saveImageAction);
        toolBarManager.add(new Separator());
    }
    
    @Override
    protected void customizeTreeViewer(TreeViewer treeViewer) {
        treeViewer.setLabelProvider(new DelegatingStyledCellLabelProvider(new GraphTreeLabelProvider()));
        // Add drop support for file
        int ops = DND.DROP_COPY | DND.DROP_MOVE;
        Transfer[] transfers = new Transfer[] {FileTransfer.getInstance()};
        treeViewer.addDropSupport(ops, transfers, new LayoutGraphViewDropAdapter(this, treeViewer));
    }
    
    @Override
    protected void setupRemainingControls(Composite parent) {
        setupGraphListViewer(parent);
        setupCanvasStack(parent);
        
        updateActionEnablement();
    }

    private void setupGraphListViewer(Composite parent) {
        graphListViewer = new TableViewer(parent);
        graphListViewer.setContentProvider(ArrayContentProvider.getInstance());
        graphListViewer.setLabelProvider(graphListLabelProvider);
        
        graphListViewer.addSelectionChangedListener(new ISelectionChangedListener() {
            @Override
            public void selectionChanged(SelectionChangedEvent event) {
                listSelectionChanged();
            }
        });
    }

    private void setupCanvasStack(Composite parent) {
        canvasStackComposite = new Composite(parent, SWT.NONE);
        
        // Setup textual graph view
        textualGraphCanvas = new Text(canvasStackComposite, SWT.H_SCROLL | SWT.V_SCROLL);
        textualGraphCanvas.setFont(JFaceResources.getFont(JFaceResources.TEXT_FONT));
        
        // Setup visual graph view
        visualGraphCanvasComposite = new ScrolledComposite(canvasStackComposite, SWT.H_SCROLL | SWT.V_SCROLL);
        visualGraphCanvasComposite.setBackground(parent.getDisplay().getSystemColor(SWT.COLOR_LIST_BACKGROUND));
        
        visualGraphCanvas = new GraphRenderingCanvas(visualGraphCanvasComposite);
        visualGraphCanvasComposite.setContent(visualGraphCanvas);
        
        // The layout can later be used to switch between the graph viewers
        canvasStackLayout = new StackLayout();
        canvasStackComposite.setLayout(canvasStackLayout);
        canvasStackLayout.topControl = null;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Events
    
    /**
     * Enables or disables our actions as required.
     */
    private void updateActionEnablement() {
        reloadFromFileAction.updateEnablement();
        saveImageAction.updateEnablement();
    }
    
}
