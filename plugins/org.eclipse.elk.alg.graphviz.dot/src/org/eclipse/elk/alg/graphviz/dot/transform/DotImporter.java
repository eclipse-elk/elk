/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.transform;

import java.util.EnumSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import org.eclipse.elk.alg.graphviz.dot.dot.Attribute;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeType;
import org.eclipse.elk.alg.graphviz.dot.dot.DotFactory;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget;
import org.eclipse.elk.alg.graphviz.dot.dot.Graph;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel;
import org.eclipse.elk.alg.graphviz.dot.dot.Node;
import org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.Statement;
import org.eclipse.elk.alg.graphviz.dot.dot.Subgraph;
import org.eclipse.elk.alg.graphviz.dot.dot.util.DotSwitch;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.MapPropertyHolder;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.util.EcoreUtil;

import com.google.common.collect.Maps;

/**
 * A transformer for Graphviz Dot.
 *
 * @author msp
 */
public class DotImporter {
    
    /** map of Graphviz node identifiers to their KNode instances. */
    private static final IProperty<Map<String, ElkNode>> NODE_ID_MAP = new Property<>("nodeIdMap");
    /** map of Graphviz port identifiers to their KPort instances. */
    private static final IProperty<Map<Pair<ElkNode, String>, ElkPort>> PORT_ID_MAP = new Property<>("portIdMap");
    /** original Graphviz statements attached to graph elements. */
    private static final IProperty<Statement> PROP_STATEMENT = new Property<>("dotTransformer.statement");
    /** original Graphviz identifiers attached to graph elements. */
    private static final IProperty<String> PROP_ID = new Property<>("dotTransformer.name");
    /** original Graphviz graph attached to parent nodes. */
    private static final IProperty<Graph> PROP_GRAPH = new Property<>("dotTransformer.graph");
    /** default node width to apply for all nodes. */
    private static final IProperty<Float> PROP_DEF_WIDTH = new Property<>("dotTransformer.defWidth");
    /** default node height to apply for all nodes. */
    private static final IProperty<Float> PROP_DEF_HEIGHT = new Property<>("dotTransformer.defHeight");
    
    /**
     * Transform the GraphViz model into a KGraph.
     * 
     * @param transData
     *            the transformation data instance that holds the source graph and is enriched with
     *            the new target graphs
     */
    public void transform(final IDotTransformationData<GraphvizModel, ElkNode> transData) {
        for (Graph graph : transData.getSourceGraph().getGraphs()) {
            ElkNode parent = ElkGraphUtil.createGraph();
            Map<String, ElkNode> nodeIdMap = Maps.newHashMap();
            transData.setProperty(NODE_ID_MAP, nodeIdMap);
            Map<Pair<ElkNode, String>, ElkPort> portIdMap = Maps.newHashMap();
            transData.setProperty(PORT_ID_MAP, portIdMap);
            transform(graph.getStatements(), parent, transData, new MapPropertyHolder(), new MapPropertyHolder());
            transData.getTargetGraphs().add(parent);
            parent.setProperty(PROP_GRAPH, graph);
        }
    }

    /**
     * Apply the layout of the target KGraphs to the original GraphViz model. This may only be used
     * on target graphs that were created by the same transformation class.
     * 
     * @param transData
     *            the transformation data instance
     */
    public void transferLayout(final IDotTransformationData<GraphvizModel, ElkNode> transData) {
        for (ElkNode layoutNode : transData.getTargetGraphs()) {
            applyLayout(layoutNode, new KVector(), layoutNode.getProperty(PROP_GRAPH));
        }
    }
    
    
    /*---------- Transformation Dot to KGraph ----------*/

    /**
     * Transform a Dot graph to a KNode.
     * 
     * @param statements
     *            a list of Dot statements
     * @param parent
     *            a KNode
     * @param transData
     *            transformation data instance
     * @param nodeProps
     *            properties that are applied to all nodes
     * @param edgeProps
     *            properties that are applied to all edges
     */
    private void transform(final List<Statement> statements, final ElkNode parent,
            final IDotTransformationData<GraphvizModel, ElkNode> transData,
            final IPropertyHolder nodeProps, final IPropertyHolder edgeProps) {
        
        DotSwitch<Object> statementSwitch = new DotSwitch<Object>() {

            public Object caseNodeStatement(final NodeStatement statement) {
                transformNode(statement, parent, transData, nodeProps);
                return null;
            }

            public Object caseEdgeStatement(final EdgeStatement statement) {
                transformEdge(statement, parent, transData, edgeProps);
                return null;
            }

            public Object caseSubgraph(final Subgraph subgraph) {
                ElkNode subElkNode = parent;
                if (subgraph.getName() != null && subgraph.getName().startsWith("cluster")) {
                    subElkNode = transformNode(subgraph.getName(), parent, transData);
                    if (subElkNode.getProperty(PROP_STATEMENT) != null) {
                        transData.log("Discarding cluster subgraph \"" + subgraph.getName()
                                + "\" since its id is already used.");
                        return null;
                    } else {
                        // the subgraph inherits all settings of its parent
                        subElkNode.copyProperties(parent);
                        subElkNode.setProperty(PROP_STATEMENT, subgraph);
                    }
                }
                MapPropertyHolder subNodeProps = new MapPropertyHolder();
                subNodeProps.copyProperties(nodeProps);
                MapPropertyHolder subEdgeProps = new MapPropertyHolder();
                subEdgeProps.copyProperties(edgeProps);
                transform(subgraph.getStatements(), subElkNode, transData, subNodeProps, subEdgeProps);
                return null;
            }

            public Object caseAttributeStatement(final AttributeStatement statement) {
                switch (statement.getType()) {
                case GRAPH:
                    for (Attribute attr : statement.getAttributes()) {
                        caseAttribute(attr);
                    }
                    break;
                case NODE:
                    for (Attribute attr : statement.getAttributes()) {
                        transformAttribute(nodeProps, attr, transData);
                    }
                    break;
                case EDGE:
                    for (Attribute attr : statement.getAttributes()) {
                        transformAttribute(edgeProps, attr, transData);
                    }
                    break;
                }
                return null;
            }

            public Object caseAttribute(final Attribute attribute) {
                if (Attributes.MARGIN.equals(attribute.getName())) {
                    ElkPadding padding = parent.getProperty(CoreOptions.PADDING);
                    if (attribute.getValue().indexOf(',') >= 0) {
                        KVector value = new KVector();
                        try {
                            value.parse(attribute.getValue());
                            padding.setLeft((float) value.x);
                            padding.setRight((float) value.x);
                            padding.setTop((float) value.y);
                            padding.setBottom((float) value.y);
                        } catch (IllegalArgumentException exception) {
                            transData.log("Discarding attribute \"" + attribute.getName()
                                    + "\" since its value could not be parsed correctly.");
                        }
                    } else {
                        try {
                            float value = Float.parseFloat(trimValue(attribute));
                            padding.setLeft(value);
                            padding.setRight(value);
                            padding.setTop(value);
                            padding.setBottom(value);
                        } catch (NumberFormatException exception) {
                            transData.log("Discarding attribute \"" + attribute.getName()
                                    + "\" since its value could not be parsed correctly.");
                        }
                    }
                } else {
                    transformAttribute(parent, attribute, transData);
                }
                return null;
            }

        };
        for (Statement statement : statements) {
            statementSwitch.doSwitch(statement);
        }
    }
    
    /**
     * Apply the property given in the attribute to a property holder.
     * 
     * @param target the property holder that receives the new property
     * @param attribute a Graphviz attribute that contains the property
     * @param transData transformation data
     */
    private void transformAttribute(final IPropertyHolder target, final Attribute attribute,
            final IDotTransformationData<GraphvizModel, ElkNode> transData) {
        
        String name = attribute.getName();
        String value = trimValue(attribute);
        try {
            if (Attributes.LAYOUT.equals(name)) {
                Command command = Command.parse(value);
                if (command != Command.INVALID) {
                    target.setProperty(CoreOptions.ALGORITHM,
                            "org.eclipse.elk.algorithm.graphviz." + command);
                } else {
                    target.setProperty(CoreOptions.ALGORITHM, value);
                }
            } else if (Attributes.WIDTH.equals(name)) {
                target.setProperty(PROP_DEF_WIDTH, Float.valueOf(value)
                        * DotExporter.DPI);
            } else if (Attributes.HEIGHT.equals(name)) {
                target.setProperty(PROP_DEF_HEIGHT, Float.valueOf(value)
                        * DotExporter.DPI);
            } else if (Attributes.FIXEDSIZE.equals(name)) {
                Boolean fixedSize = Boolean.valueOf(value);
                target.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS,
                        fixedSize ? SizeConstraint.fixed() : EnumSet.of(SizeConstraint.MINIMUM_SIZE));
            } else if (Attributes.NODESEP.equals(name)) {
                target.setProperty(CoreOptions.SPACING_NODE_NODE, Double.valueOf(value));
            } else if (Attributes.PACK.equals(name)) {
                target.setProperty(CoreOptions.SEPARATE_CONNECTED_COMPONENTS, Boolean.valueOf(value));
            } else if (Attributes.PAD.equals(name)) {
                if (value.indexOf(',') >= 0) {
                    KVector pad = new KVector();
                    pad.parse(value);
                    ElkPadding padding = new ElkPadding(pad.x / 2, pad.y / 2);
                    target.setProperty(CoreOptions.PADDING, padding);
                } else {
                    target.setProperty(CoreOptions.PADDING, new ElkPadding(Float.valueOf(value)));
                }
            } else if (Attributes.RANKDIR.equals(name)) {
                if (value.equals("TB")) {
                    target.setProperty(CoreOptions.DIRECTION, Direction.DOWN);
                } else if (value.equals("BT")) {
                    target.setProperty(CoreOptions.DIRECTION, Direction.UP);
                } else if (value.equals("LR")) {
                    target.setProperty(CoreOptions.DIRECTION, Direction.RIGHT);
                } else if (value.equals("RL")) {
                    target.setProperty(CoreOptions.DIRECTION, Direction.LEFT);
                }
            } else if (Attributes.SPLINES.equals(name)) {
                if (value.equals("spline") || value.equals("true")) {
                    target.setProperty(CoreOptions.EDGE_ROUTING, EdgeRouting.SPLINES);
                } else if (value.equals("polyline") || value.equals("line")
                        || value.equals("false")) {
                    target.setProperty(CoreOptions.EDGE_ROUTING, EdgeRouting.POLYLINE);
                } else if (value.equals("ortho")) {
                    target.setProperty(CoreOptions.EDGE_ROUTING, EdgeRouting.ORTHOGONAL);
                }
            } else if (Attributes.START.equals(name)) {
                if (value.startsWith("random")) {
                    int random = 0;
                    try {
                        random = Integer.parseInt(value.substring("random".length()));
                    } catch (NumberFormatException e) {
                        // ignore exception
                    }
                    target.setProperty(CoreOptions.RANDOM_SEED, random);
                }
            } else if (Attributes.WEIGHT.equals(name)) {
                target.setProperty(CoreOptions.PRIORITY, (int) Float.parseFloat(value));
            }
        } catch (NumberFormatException exception) {
            transData.log("Discarding attribute " + attribute.getName() + "="
                    + attribute.getValue() + " since its value could not be parsed correctly.");
        } catch (IllegalArgumentException exception) {
            transData.log("Discarding attribute " + attribute.getName() + "="
                    + attribute.getValue() + " since its value could not be parsed correctly.");
        }
    }
    
    /**
     * Transforms a node.
     * 
     * @param statement a node statement
     * @param parent the parent node
     * @param transData the transformation data instance
     * @param defaultProps default values for node options
     */
    private void transformNode(final NodeStatement statement, final ElkNode parent,
            final IDotTransformationData<GraphvizModel, ElkNode> transData,
            final IPropertyHolder defaultProps) {
        
        ElkNode elknode = transformNode(statement.getNode().getName(), parent, transData);
        if (elknode.getProperty(PROP_STATEMENT) != null) {
            transData.log("Discarding node \"" + statement.getNode().getName()
                    + "\" since its id is already used.");
        } else {
            elknode.copyProperties(defaultProps);
            elknode.setProperty(PROP_STATEMENT, statement);
            Float defWidth = defaultProps.getProperty(PROP_DEF_WIDTH);
            if (defWidth != null) {
                elknode.setWidth(defWidth);
            }
            Float defHeight = defaultProps.getProperty(PROP_DEF_HEIGHT);
            if (defHeight != null) {
                elknode.setHeight(defHeight);
            }
            
            // evaluate attributes for the new node
            for (Attribute attr : statement.getAttributes()) {
                String value = trimValue(attr);
                try {
                    if (Attributes.LABEL.equals(attr.getName())) {
                        ElkGraphUtil.createLabel(value, elknode);
                    } else if (Attributes.POS.equals(attr.getName())) {
                        KVector pos = new KVector();
                        pos.parse(value);
                        pos.scale(DotExporter.DPI);
                        elknode.setLocation(pos.x, pos.y);
                    } else if (Attributes.WIDTH.equals(attr.getName())) {
                        elknode.setWidth(Float.parseFloat(value) * DotExporter.DPI);
                    } else if (Attributes.HEIGHT.equals(attr.getName())) {
                        elknode.setHeight(Float.parseFloat(value) * DotExporter.DPI);
                    } else {
                        transformAttribute(elknode, attr, transData);
                    }
                } catch (NumberFormatException exception) {
                    transData.log("Discarding attribute \"" + attr.getName()
                            + "\" for node \"" + statement.getNode().getName()
                            + "\" since its value could not be parsed correctly.");
                } catch (IllegalArgumentException exception) {
                    transData.log("Discarding attribute \"" + attr.getName()
                            + "\" for node \"" + statement.getNode().getName()
                            + "\" since its value could not be parsed correctly.");
                }
            }
        }
    }
    
    /** the default width and height value for nodes. */
    private static final float DEF_WIDTH = 10.0f;
    
    /**
     * Transforms a single node, if not already done before.
     * 
     * @param nodeId a node identifier
     * @param parent the parent where the new KNode is stored
     * @param transData the transformation data instance
     * @return a KNode instance
     */
    private ElkNode transformNode(final String nodeId, final ElkNode parent,
            final IDotTransformationData<GraphvizModel, ElkNode> transData) {
        
        Map<String, ElkNode> nodeIdMap = transData.getProperty(NODE_ID_MAP);
        ElkNode elknode = nodeIdMap.get(nodeId);
        if (elknode == null) {
            elknode = ElkGraphUtil.createNode(parent);
            elknode.setWidth(DEF_WIDTH);
            elknode.setHeight(DEF_WIDTH);
            if (nodeId != null) {
                nodeIdMap.put(nodeId, elknode);
                elknode.setProperty(PROP_ID, nodeId);
            }
        }
        return elknode;
    }
    
    /**
     * Transforms a single port, if not already done before.
     * 
     * @param portId a port identifier
     * @param node the node to which the new KPort belongs
     * @param transData the transformation data instance
     * @return a KPort instance
     */
    private ElkPort transformPort(final String portId, final ElkNode node,
            final IDotTransformationData<GraphvizModel, ElkNode> transData) {
        
        Map<Pair<ElkNode, String>, ElkPort> portIdMap = transData.getProperty(PORT_ID_MAP);
        Pair<ElkNode, String> key = new Pair<>(node, portId);
        ElkPort elkport = portIdMap.get(key);
        if (elkport == null) {
            elkport = ElkGraphUtil.createPort(node);
            if (portId != null) {
                portIdMap.put(key, elkport);
            }
        }
        return elkport;
    }
    
    /**
     * Transforms an edge.
     * 
     * @param statement an edge statement
     * @param parent the parent node
     * @param transData the transformation data instance
     * @param defaultProps default values for edge options
     */
    private void transformEdge(final EdgeStatement statement, final ElkNode parent,
            final IDotTransformationData<GraphvizModel, ElkNode> transData,
            final IPropertyHolder defaultProps) {
        
        String sourceName = statement.getSourceNode().getName();
        ElkNode source = transformNode(sourceName, parent, transData);
        ElkPort sourcePort = null;
        if (statement.getSourceNode().getPort() != null) {
            String portName = statement.getSourceNode().getPort().getCompass_pt();
            if (portName == null) {
                portName = statement.getSourceNode().getPort().getName();
            }
            sourcePort = transformPort(portName, source, transData);
        }
        ListIterator<EdgeTarget> targetIter = statement.getEdgeTargets().listIterator();
        while (targetIter.hasNext()) {
            EdgeTarget edgeTarget = targetIter.next();
            ElkEdge elkedge = ElkGraphUtil.createEdge(null);
            if (sourcePort != null) {
                elkedge.getSources().add(sourcePort);
            } else {
                elkedge.getSources().add(source);
            }
            ElkNode target;
            ElkPort targetPort = null;
            Node edgeTargetNode = edgeTarget.getTargetnode();
            String targetName;
            if (edgeTargetNode == null) {
                targetName = edgeTarget.getTargetSubgraph().getName();
                target = transformNode(targetName, parent, transData);
            } else {
                targetName = edgeTargetNode.getName();
                target = transformNode(targetName, parent, transData);
                if (edgeTargetNode.getPort() != null) {
                    String portName = edgeTargetNode.getPort().getCompass_pt();
                    if (portName == null) {
                        portName = edgeTargetNode.getPort().getName();
                    }
                    targetPort = transformPort(portName, target, transData);
                }
            }
            
            if (targetPort != null) {
                elkedge.getTargets().add(targetPort);
            } else {
                elkedge.getTargets().add(target);
            }
            
            elkedge.copyProperties(defaultProps);
            if (targetIter.previousIndex() == 0) {
                // this is the first target - just store the edge statement
                elkedge.setProperty(PROP_STATEMENT, statement);
            } else {
                // the edge statement has more that one target - create a copy
                EdgeStatement newStatement = DotFactory.eINSTANCE.createEdgeStatement();
                Node sourceNode = DotFactory.eINSTANCE.createNode();
                sourceNode.setName(sourceName);
                newStatement.setSourceNode(sourceNode);
                targetIter.remove();
                newStatement.getEdgeTargets().add(edgeTarget);
                for (Attribute attr : statement.getAttributes()) {
                    newStatement.getAttributes().add(EcoreUtil.copy(attr));
                }
                elkedge.setProperty(PROP_STATEMENT, newStatement);
            }
            
            // evaluate attributes for the new edge
            for (Attribute attr : statement.getAttributes()) {
                String value = trimValue(attr);
                if (Attributes.LABEL.equals(attr.getName())) {
                    ElkLabel label = ElkGraphUtil.createLabel(value, elkedge);
                    label.setProperty(CoreOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.CENTER);
                } else if (Attributes.HEADLABEL.equals(attr.getName())) {
                    ElkLabel label = ElkGraphUtil.createLabel(value, elkedge);
                    label.setProperty(CoreOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.HEAD);
                } else if (Attributes.TAILLABEL.equals(attr.getName())) {
                    ElkLabel label = ElkGraphUtil.createLabel(value, elkedge);
                    label.setProperty(CoreOptions.EDGE_LABELS_PLACEMENT, EdgeLabelPlacement.TAIL);
                } else {
                    transformAttribute(elkedge, attr, transData);
                }
            }
            
            // the edge target is the source for the next edge target
            source = target;
            sourceName = targetName;
            sourcePort = targetPort;
        }
    }
    

    /*---------- Layout Transfer KGraph to Dot ----------*/

    /**
     * Apply layout to a parent node and its children.
     * 
     * @param parent a parent node
     * @param offset the node's offset in the graph
     * @param graph the Graphviz graph
     */
    private void applyLayout(final ElkNode parent, final KVector offset, final Graph graph) {
        for (ElkNode elknode : parent.getChildren()) {
            Statement statement = elknode.getProperty(PROP_STATEMENT);
            if (statement == null) {
                // the node was only declared implicitly - create an explicit declaration
                NodeStatement stm = DotFactory.eINSTANCE.createNodeStatement();
                Node node = DotFactory.eINSTANCE.createNode();
                node.setName(elknode.getProperty(PROP_ID));
                stm.setNode(node);
                graph.getStatements().add(stm);
                statement = stm;
            }
            if (statement instanceof NodeStatement) {
                List<Attribute> attributes = ((NodeStatement) statement).getAttributes();
                // transfer node position
                removeAttributes(attributes, Attributes.POS);
                double xpos = elknode.getX() + elknode.getWidth() / 2 + offset.x;
                double ypos = elknode.getY() + elknode.getHeight() / 2 + offset.y;
                String posString = "\"" + Double.toString(xpos) + "," + Double.toString(ypos) + "\"";
                attributes.add(DotExporter.createAttribute(Attributes.POS, posString));
                // transfer node size
                removeAttributes(attributes, Attributes.WIDTH);
                attributes.add(DotExporter.createAttribute(Attributes.WIDTH, elknode.getWidth() / DotExporter.DPI));
                removeAttributes(attributes, Attributes.HEIGHT);
                attributes.add(DotExporter.createAttribute(Attributes.HEIGHT, elknode.getHeight() / DotExporter.DPI));
            } else if (statement instanceof Subgraph) {
                applyLayout(elknode, new KVector(offset).add(elknode.getX(), elknode.getY()), graph);
            }
            
            for (ElkEdge elkedge : ElkGraphUtil.allOutgoingEdges(elknode)) {
                applyLayout(elkedge, offset, graph);
            }
        }
        
        // transfer graph size to bounding box
        List<Statement> statements;
        Statement graphStm = parent.getProperty(PROP_STATEMENT);
        if (graphStm instanceof Subgraph) {
            statements = ((Subgraph) graphStm).getStatements();
        } else {
            statements = graph.getStatements();
        }
        removeGraphAttributes(statements, Attributes.BOUNDINGBOX);
        String bbString = "\"0,0," + Double.toString(parent.getWidth())
                + "," + Double.toString(parent.getHeight()) + "\"";
        statements.add(DotExporter.createAttribute(
                Attributes.BOUNDINGBOX, bbString));
    }
    
    /**
     * Apply layout to an edge and its labels.
     * 
     * @param edge an edge
     * @param offset its offset in the graph
     * @param graph the Graphviz graph
     */
    private void applyLayout(final ElkEdge edge, final KVector offset, final Graph graph) {
        EdgeStatement edgeStatement = (EdgeStatement) edge.getProperty(PROP_STATEMENT);
        if (edgeStatement.eContainer() == null) {
            // this can happen when an edge with multiple target declarations was found
            graph.getStatements().add(edgeStatement);
        }
        
        // transfer edge bend points and source / target points
        List<Attribute> attributes = edgeStatement.getAttributes();
        removeAttributes(attributes, Attributes.POS);
        
        if (!edge.getSections().isEmpty()) {
            StringBuilder bendpointString = new StringBuilder("\"");
            KVectorChain vectorChain = ElkUtil.createVectorChain(edge.getSections().get(0));
            ListIterator<KVector> chainIter = vectorChain.listIterator();
            while (chainIter.hasNext()) {
                KVector point = chainIter.next().add(offset);
                bendpointString.append(point.x);
                bendpointString.append(',');
                bendpointString.append(point.y);
                if (chainIter.hasNext()) {
                    bendpointString.append(' ');
                }
            }
            bendpointString.append('\"');
            attributes.add(DotExporter.createAttribute(Attributes.POS,
                    bendpointString.toString()));
        }
        
        // transfer label positions
        for (ElkLabel label : edge.getLabels()) {
            String attrKey = null;
            switch (label.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT)) {
            case CENTER:
                attrKey = Attributes.LABELPOS;
                break;
            case HEAD:
                attrKey = Attributes.HEADLP;
                break;
            case TAIL:
                attrKey = Attributes.TAILLP;
                break;
            }
            if (attrKey != null) {
                removeAttributes(attributes, attrKey);
                double xpos = label.getX() + label.getWidth() / 2 + offset.x;
                double ypos = label.getY() + label.getHeight() / 2 + offset.y;
                String posString = "\"" + Double.toString(xpos)
                        + "," + Double.toString(ypos) + "\"";
                attributes.add(DotExporter.createAttribute(attrKey, posString));
            }
        }
    }
    
    /**
     * Removes all graph-related occurrences of the given key from a statement list.
     * 
     * @param statements a list of statements
     * @param key an attribute key
     */
    private void removeGraphAttributes(final List<Statement> statements, final String key) {
        ListIterator<Statement> stmIter = statements.listIterator();
        while (stmIter.hasNext()) {
            Statement stm = stmIter.next();
            if (stm instanceof Attribute) {
                if (key.equals(((Attribute) stm).getName())) {
                    stmIter.remove();
                }
            } else if (stm instanceof AttributeStatement) {
                AttributeStatement attrStatement = (AttributeStatement) stm;
                if (attrStatement.getType() == AttributeType.GRAPH) {
                    removeAttributes(attrStatement.getAttributes(), key);
                    if (attrStatement.getAttributes().isEmpty()) {
                        stmIter.remove();
                    }
                }
            }
        }
    }
    
    /**
     * Removes all occurrences of the given key from an attributes list.
     * 
     * @param attributes a list of attributes
     * @param key an attribute key
     */
    private void removeAttributes(final List<Attribute> attributes, final String key) {
        ListIterator<Attribute> attrIter = attributes.listIterator();
        while (attrIter.hasNext()) {
            Attribute attr = attrIter.next();
            if (key.equals(attr.getName())) {
                attrIter.remove();
            }
        }
    }
    
    /**
     * Remove the quote characters leading and trailing the attribute value.
     * 
     * @param attribute an attribute
     * @return a trimmed value of the attribute
     */
    private String trimValue(final Attribute attribute) {
        String value = attribute.getValue();
        if (value == null) {
            value = "true";
        } else if (value.startsWith("\"") && value.endsWith("\"")) {
            return value.substring(1, value.length() - 1).trim();
        }
        return value;
    }
    
}
