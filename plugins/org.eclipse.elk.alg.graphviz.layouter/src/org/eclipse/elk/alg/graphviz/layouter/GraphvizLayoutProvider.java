/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.layouter;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Map;

import org.eclipse.core.runtime.Platform;
import org.eclipse.elk.alg.graphviz.dot.GraphvizDotStandaloneSetup;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel;
import org.eclipse.elk.alg.graphviz.dot.transform.Command;
import org.eclipse.elk.alg.graphviz.dot.transform.DotExporter;
import org.eclipse.elk.alg.graphviz.dot.transform.DotResourceSetProvider;
import org.eclipse.elk.alg.graphviz.dot.transform.DotTransformationData;
import org.eclipse.elk.alg.graphviz.dot.transform.IDotTransformationData;
import org.eclipse.elk.alg.graphviz.layouter.GraphvizTool.Cleanup;
import org.eclipse.elk.alg.graphviz.layouter.preferences.GraphvizLayouterPreferenceStoreAccess;
import org.eclipse.elk.alg.graphviz.layouter.util.ForkedOutputStream;
import org.eclipse.elk.alg.graphviz.layouter.util.ForwardingInputStream;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.WrappedException;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.common.util.URI;
import org.eclipse.emf.ecore.resource.Resource.Diagnostic;
import org.eclipse.xtext.resource.SaveOptions;
import org.eclipse.xtext.resource.XtextResource;
import org.eclipse.xtext.resource.XtextResourceSet;

import com.google.inject.Injector;

/**
 * Layout provider for the Graphviz layout tool.
 * The actual Graphviz layout that is applied is determined by the parameter
 * passed in the {@link #initialize(String)} method.
 * 
 * @author msp
 */
public class GraphvizLayoutProvider extends AbstractLayoutProvider {

    /** preference constant for determining whether to reuse a single Graphviz process. */
    public static final String PREF_GRAPHVIZ_REUSE_PROCESS = "graphviz.reuseProcess";
    /** default setting of above defined preference. */
    public static final boolean REUSE_PROCESS_DEFAULT = true;

    /** the serial call number for usage in debug mode. */
    private static int serialCallNo = 0;
    
    /** command passed to the layouter. */
    private Command command = Command.INVALID;
    /** the Graphviz process pool. */
    private GraphvizTool graphvizTool;
    /** the Graphviz Dot format handler. */
    private DotResourceSetProvider dotResourceSetProvider;
    /** the call number for the current execution. */
    private int myCallNo;
    /** lazily created injector for creating required resource set creators. */
    private Injector injector;

    @Override
    public void initialize(final String parameter) {
        command = Command.valueOf(parameter);
        graphvizTool = new GraphvizTool(command);
        
        // the dot format handler is indirectly fetched in order to ensure proper injection (if we're
        // inside Eclipse, use the GraphFormatsService to retrieve the handler; otherwise, use an
        // injector to retrieve an instance)
        dotResourceSetProvider = getInjector().getInstance(DotResourceSetProvider.class);
    }

    @Override
    public void dispose() {
        graphvizTool.cleanup(Cleanup.STOP);
    }
    
    /**
     * Returns the injector, creating a new instance if none was created yet.
     * 
     * @return the injector.
     */
    private Injector getInjector() {
        if (injector == null) {
            injector = new GraphvizDotStandaloneSetup().createInjectorAndDoEMFRegistration();
        }
        
        return injector;
    }

    @Override
    public void layout(final ElkNode parentNode, final IElkProgressMonitor progressMonitor) {
        if (command == Command.INVALID) {
            throw new IllegalStateException("The Graphviz layout provider is not initialized.");
        }
        progressMonitor.begin("Graphviz layout (" + command + ")", 2);
        if (parentNode.getChildren().isEmpty()) {
            // return if there is nothing in this node
            progressMonitor.done();
            return;
        }
        boolean debugMode = parentNode.getProperty(CoreOptions.DEBUG_MODE);
        myCallNo = ++serialCallNo;

        // start the graphviz process, or retrieve the previously used process
        graphvizTool.initialize();

        // create an Xtext resource set for parsing and serialization
        XtextResourceSet resourceSet = (XtextResourceSet) dotResourceSetProvider.createResourceSet();
        
        // create the dot exporter we'll be using
        DotExporter dotExporter = new LayoutDotExporter();

        // translate the KGraph to Graphviz and write to the process
        IDotTransformationData<ElkNode, GraphvizModel> transData = new DotTransformationData<ElkNode, GraphvizModel>();
        transData.setSourceGraph(parentNode);
        transData.setProperty(DotExporter.COMMAND, command);
        
        dotExporter.transform(transData);
        GraphvizModel graphvizInput = transData.getTargetGraphs().get(0);
        writeDotGraph(graphvizInput, progressMonitor.subTask(1), debugMode, resourceSet);

        try {
            // read Graphviz output and apply layout information to the KGraph
            GraphvizModel graphvizOutput = readDotGraph(progressMonitor.subTask(1),
                    debugMode, resourceSet);
            transData.getTargetGraphs().set(0, graphvizOutput);
            dotExporter.transferLayout(transData);
        } finally {
            boolean reuseProcess = false;
            if (Platform.isRunning()) {
                reuseProcess = GraphvizLayouterPreferenceStoreAccess.getUISaveBoolean(
                        PREF_GRAPHVIZ_REUSE_PROCESS, REUSE_PROCESS_DEFAULT);
            }
            graphvizTool.cleanup(reuseProcess ? Cleanup.NORMAL : Cleanup.STOP);
            progressMonitor.done();
        }
    }

    /**
     * Writes a serialized version of the Graphviz model to the given output stream.
     * 
     * @param graphvizModel
     *            Graphviz model to serialize
     * @param monitor
     *            a monitor to which progress is reported
     * @param debugMode
     *            whether debug mode is active
     * @param resourceSet
     *            the resource set for serialization
     */
    private void writeDotGraph(final GraphvizModel graphvizModel,
            final IElkProgressMonitor monitor, final boolean debugMode,
            final XtextResourceSet resourceSet) {
        monitor.begin("Serialize model", 1);
        OutputStream outputStream = graphvizTool.input();
        // enable debug output if needed
        FileOutputStream debugStream = null;
        if (debugMode) {
            try {
                String path = System.getProperty("user.home");
                if (path.endsWith(File.separator)) {
                    path += "tmp" + File.separator + "graphviz";
                } else {
                    path += File.separator + "tmp" + File.separator + "graphviz";
                }
                new File(path).mkdirs();
                debugStream = new FileOutputStream(new File(path + File.separator
                        + debugFileBase() + "-in.dot"));
                outputStream = new ForkedOutputStream(outputStream, debugStream);
            } catch (Exception exception) {
                System.out.println("GraphvizLayouter: Could not initialize debug output: "
                        + exception.getMessage());
            }
        }

        try {
            XtextResource resource = (XtextResource) resourceSet.createResource(
                    URI.createURI("output.graphviz_dot"));
            resource.getContents().add(graphvizModel);
            
            /* KIPRA-1498
             * We disable formatting and validation when saving the resource. Enabling it lead to
             * possible ConcurrentModificationExceptions in Xtext.
             */
            Map<Object, Object> saveOptions =
                    SaveOptions.newBuilder().noValidation().getOptions().toOptionsMap();
            resource.save(outputStream, saveOptions);
            
            outputStream.write('\n');
            outputStream.flush();
        } catch (IOException exception) {
            graphvizTool.cleanup(Cleanup.ERROR);
            throw new WrappedException("Failed to send the graph to Graphviz.", exception);
        } finally {
            if (debugStream != null) {
                try {
                    debugStream.close();
                } catch (IOException exception) {
                    // ignore exception
                }
            }
            monitor.done();
        }
    }

    /**
     * Reads and parses a serialized Graphviz model.
     * 
     * @param monitor
     *            a monitor to which progress is reported
     * @param debugMode
     *            whether debug mode is active
     * @param resourceSet
     *            the resoure set for parsing
     * @return an instance of the parsed graphviz model
     */
    private GraphvizModel readDotGraph(final IElkProgressMonitor monitor,
            final boolean debugMode, final XtextResourceSet resourceSet) {
        monitor.begin("Parse output", 1);
        InputStream inputStream = graphvizTool.output();
        // enable debug output if needed
        FileOutputStream debugStream = null;
        if (debugMode) {
            try {
                String path = System.getProperty("user.home");
                if (path.endsWith(File.separator)) {
                    path += "tmp" + File.separator + "graphviz";
                } else {
                    path += File.separator + "tmp" + File.separator + "graphviz";
                }
                new File(path).mkdirs();
                debugStream = new FileOutputStream(new File(path + File.separator
                        + debugFileBase() + "-out.dot"));
                inputStream = new ForwardingInputStream(inputStream, debugStream);
            } catch (Exception exception) {
                System.out.println("GraphvizLayouter: Could not initialize debug output: "
                        + exception.getMessage());
            }
        }

        // parse the output stream of the dot process
        XtextResource resource = (XtextResource) resourceSet.createResource(
                URI.createURI("input.graphviz_dot"));
        try {
            resource.load(inputStream, null);
        } catch (IOException exception) {
            graphvizTool.cleanup(Cleanup.ERROR);
            throw new WrappedException("Failed to read Graphviz output.", exception);
        } finally {
            if (debugStream != null) {
                try {
                    debugStream.close();
                } catch (IOException exception) {
                    // ignore exception
                }
            }
        }

        // analyze errors and retrieve parse result
        if (!resource.getErrors().isEmpty()) {
            StringBuilder errorString = new StringBuilder("Errors in Graphviz output:");
            for (Diagnostic diagnostic : resource.getErrors()) {
                errorString.append("\n" + diagnostic.getLine() + ": " + diagnostic.getMessage());
            }
            graphvizTool.cleanup(Cleanup.ERROR);
            throw new GraphvizException(errorString.toString());
        }
        GraphvizModel graphvizModel = (GraphvizModel) resource.getParseResult().getRootASTElement();
        if (graphvizModel == null || graphvizModel.getGraphs().isEmpty()) {
            graphvizTool.cleanup(Cleanup.ERROR);
            throw new GraphvizException("No output from the Graphviz process."
                    + " Try increasing the timeout value in the Eclipse Diagram Layout preferences.");
        }

        monitor.done();
        return graphvizModel;
    }
    
    /**
     * Return the base name for debug files.
     * 
     * @return the base name for debug files
     */
    private String debugFileBase() {
        String no = Integer.toString(myCallNo);
        switch (no.length()) {
        case 1:
            return "debug00" + no;
        case 2:
            return "debug0" + no;
        default:
            return "debug" + no;
        }
    }
    
}
