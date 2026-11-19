/*******************************************************************************
 * Copyright (c) 2010, 2020 Kiel University and others.
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

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.eclipse.elk.alg.graphviz.dot.dot.Attribute;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeType;
import org.eclipse.elk.alg.graphviz.dot.dot.DotFactory;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget;
import org.eclipse.elk.alg.graphviz.dot.dot.Graph;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphType;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel;
import org.eclipse.elk.alg.graphviz.dot.dot.Node;
import org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.Statement;
import org.eclipse.elk.alg.graphviz.dot.dot.Subgraph;
import org.eclipse.elk.core.UnsupportedGraphException;
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
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

/**
 * This class implements a transformation from the KGraph metamodel to the Dot metamodel.
 * Furthermore it contains functionality to apply layout information attached to a Dot model to a
 * KGraph model.
 */
public class DotExporter {

    /** small default value for minimal spacing. */
    public static final double DEF_SPACING_SMALL = 20.0;
    /** large default value for minimal spacing. */
    public static final double DEF_SPACING_LARGE = 40.0;
    /** extra-large default value for minimal spacing. */
    public static final double DEF_SPACING_XLARGE = 60.0;
    /** dots per inch specification, needed by Graphviz for some values. */
    public static final float DPI = 72.0f;

    /** the Graphviz command to use for transformation. */
    public static final IProperty<Command> COMMAND = new Property<Command>(
            "dotExporter.command", Command.DOT);
    /** whether edge identifiers should be generated or not. */
    public static final IProperty<Boolean> USE_EDGE_IDS = new Property<Boolean>(
            "dotExporter.useEdgeIds", false);
    /** whether to transform all hierarchy levels of the graph. */
    public static final IProperty<Boolean> HIERARCHY = new Property<Boolean>(
            "dotExporter.hierarchy", false);
    /** whether to transform the layout information of nodes. */
    public static final IProperty<Boolean> TRANSFORM_NODE_LAYOUT = new Property<Boolean>(
            "dotExporter.transformNodeLayout", true);
    /** whether to transform the layout information of edges. */
    public static final IProperty<Boolean> TRANSFORM_EDGE_LAYOUT = new Property<Boolean>(
            "dotExporter.transformEdgeLayout", true);
    /** whether to transform node labels. */
    public static final IProperty<Boolean> TRANSFORM_NODE_LABELS = new Property<Boolean>(
            "dotExporter.transformNodeLabels", true);
    /** whether to adapt port positions to edge end points. */
    public static final IProperty<Boolean> ADAPT_PORT_POSITIONS = new Property<Boolean>(
            "dotExporter.adaptPortPositions", false);
    
    /** default multiplier for font sizes. */
    private static final double FONT_SIZE_MULT = 1.4;
    /** set of delimiters used to parse attribute values. */
    private static final String ATTRIBUTE_DELIM = "\", \t\n\r";

    /** maps each identifier of a graph element to the instance of the element. */
    private static final IProperty<BiMap<String, ElkGraphElement>> GRAPH_ELEMS =
            new Property<>("dotExporter.graphElemMap");
    /** indicates whether splines are used. */
    protected static final IProperty<Boolean> USE_SPLINES = new Property<>("dotExporter.useSplines", false);
    /** the next node identifier to use. */
    private static final IProperty<Integer> NEXT_NODE_ID = new Property<>("dotExporter.nextNodeId", 1);
    /** the next edge identifier to use. */
    private static final IProperty<Integer> NEXT_EDGE_ID = new Property<>("dotExporter.nextEdgeId", 1);
    /** cluster dummy node identifier attached to a node. */
    private static final IProperty<String> CLUSTER_DUMMY = new Property<>("dotExporter.clusterDummy");

    /**
     * Transforms the KGraph instance to a Dot instance using the given command.
     * 
     * @param transData the transformation data instance
     */
    public void transform(final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        BiMap<String, ElkGraphElement> graphElems = HashBiMap.create();
        transData.setProperty(GRAPH_ELEMS, graphElems);
        ElkNode elkgraph = transData.getSourceGraph();
        GraphvizModel graphvizModel = DotFactory.eINSTANCE.createGraphvizModel();
        Graph graph = DotFactory.eINSTANCE.createGraph();
        graph.setType(GraphType.DIGRAPH);
        graphvizModel.getGraphs().add(graph);
        transformNodes(elkgraph, graph.getStatements(), new KVector(), transData);
        transformEdges(elkgraph, graph.getStatements(), transData);
        transData.getTargetGraphs().add(graphvizModel);
    }

    /**
     * Applies the layout information attached to the given Dot instance to the KGraph instance
     * using the mapping created by a previous call to {@code transform}. Has to be called after a
     * call to {@code transform}.
     * 
     * @param transData the transformation data instance
     */
    public void transferLayout(final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        ElkPadding padding = transData.getSourceGraph().getProperty(CoreOptions.PADDING);
        Graph graph = transData.getTargetGraphs().get(0).getGraphs().get(0);
        
        // process nodes and subgraphs
        KVector baseOffset = new KVector();
        applyLayout(transData.getSourceGraph(), graph.getStatements(), baseOffset, padding, transData);
        
        // finally process the edges
        LinkedList<Statement> statements = new LinkedList<Statement>(graph.getStatements());
        KVector edgeOffset = baseOffset.add(padding.getLeft(), padding.getTop());
        while (!statements.isEmpty()) {
            Statement statement = statements.removeFirst();
            if (statement instanceof EdgeStatement) {
                applyEdgeLayout((EdgeStatement) statement, edgeOffset, transData);
            } else if (statement instanceof Subgraph) {
                statements.addAll(((Subgraph) statement).getStatements());
            }
        }
    }

    /**
     * Transform the child nodes of the given parent node.
     * 
     * @param parent a parent node
     * @param statements the list to which new statements are added
     * @param offset offset of the parent node in the whole graph
     * @param transData transformation data
     */
    private void transformNodes(final ElkNode parent, final List<Statement> statements,
            final KVector offset, final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        // set attributes for the whole graph
        setGraphAttributes(statements, parent, transData);

        // create nodes and subgraphs
        boolean hierarchy = transData.getProperty(HIERARCHY);
        boolean transformNodeLayout = transData.getProperty(TRANSFORM_NODE_LAYOUT);
        boolean transformNodeLabels = transData.getProperty(TRANSFORM_NODE_LABELS);
        for (ElkNode childNode : parent.getChildren()) {
            NodeStatement nodeStatement = DotFactory.eINSTANCE.createNodeStatement();
            List<Attribute> attributes = nodeStatement.getAttributes();
            String nodeID;
            // if hierarchy mode is active, create a subgraph, else a regular node
            if (hierarchy && !childNode.getChildren().isEmpty()) {
                String clusterNodeID = getNodeID(childNode, NodeType.CLUSTER, transData);
                Subgraph subgraph = DotFactory.eINSTANCE.createSubgraph();
                subgraph.setName(clusterNodeID);
                statements.add(subgraph);
                
                // transform child nodes recursively
                ElkPadding padding = childNode.getProperty(CoreOptions.PADDING);
                double subgraphx = childNode.getX() + padding.getLeft();
                double subgraphy = childNode.getY() + padding.getTop();
                transformNodes(childNode, subgraph.getStatements(),
                        new KVector(offset).add(subgraphx, subgraphy), transData);
                
                // create a dummy node for compound edges
                nodeID = getNodeID(childNode, NodeType.DUMMY, transData);
                attributes.add(createAttribute(Attributes.STYLE, "invis"));
                attributes.add(createAttribute(Attributes.WIDTH, 0));
                attributes.add(createAttribute(Attributes.HEIGHT, 0));
                subgraph.getStatements().add(nodeStatement);
            } else {
                nodeID = getNodeID(childNode, NodeType.NODE, transData);
                // set width and height
                ElkUtil.resizeNode(childNode);
                if (childNode.getWidth() > 0) {
                    attributes.add(createAttribute(Attributes.WIDTH, childNode.getWidth() / DPI));
                }
                if (childNode.getHeight() > 0) {
                    attributes.add(createAttribute(Attributes.HEIGHT, childNode.getHeight() / DPI));
                }
                if (transformNodeLabels && !childNode.getLabels().isEmpty()
                        && childNode.getLabels().get(0).getText().length() > 0) {
                    
                    attributes.add(createAttribute(Attributes.LABEL, createString(
                            childNode.getLabels().get(0).getText())));
                }
                // add node position if interactive layout is chosen
                if (transformNodeLayout && (childNode.getX() != 0 || childNode.getY() != 0)) {
                    double xpos = (childNode.getX() + childNode.getWidth() / 2 + offset.x);
                    double ypos = (childNode.getY() + childNode.getHeight() / 2 + offset.y);
                    String posString = "\"" + Double.toString(xpos)
                            + "," + Double.toString(ypos) + "\"";
                    attributes.add(createAttribute(Attributes.POS, posString));
                }
                statements.add(nodeStatement);
            }
            Node node = DotFactory.eINSTANCE.createNode();
            node.setName(nodeID);
            nodeStatement.setNode(node);
        }
    }

    /**
     * Transform the edges of the given parent node.
     * 
     * @param parent a parent node
     * @param statements the list to which new statements are added
     * @param transData transformation data
     */
    private void transformEdges(final ElkNode parent, final List<Statement> statements,
            final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        boolean hierarchy = transData.getProperty(HIERARCHY);
        boolean transformEdgeLayout = transData.getProperty(TRANSFORM_EDGE_LAYOUT);
        Direction direction = parent.getProperty(CoreOptions.DIRECTION);
        boolean vertical = direction == Direction.DOWN || direction == Direction.UP
                || direction == Direction.UNDEFINED;
        LinkedList<ElkNode> nodes = new LinkedList<>(parent.getChildren());
        BiMap<ElkGraphElement, String> nodeIds = transData.getProperty(GRAPH_ELEMS).inverse();
        
        while (!nodes.isEmpty()) {
            ElkNode source = nodes.removeFirst();
            for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(source)) {
                // We don't support hyperedges
                if (edge.isHyperedge()) {
                    throw new UnsupportedGraphException("Hyperedges are not supported.");
                }
                
                ElkNode target = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0));
                
                // cross-hierarchy edges are considered only if hierarchy mode is active
                if (source.getParent() == target.getParent()
                        || hierarchy && isInsideGraph(target, transData.getSourceGraph())) {
                    
                    EdgeStatement edgeStatement = DotFactory.eINSTANCE.createEdgeStatement();
                    List<Attribute> attributes = edgeStatement.getAttributes();
                    
                    // set source node or cluster
                    Node sourceNode = DotFactory.eINSTANCE.createNode();
                    if (hierarchy && !source.getChildren().isEmpty()) {
                        sourceNode.setName(source.getProperty(CLUSTER_DUMMY));
                        attributes.add(createAttribute(Attributes.LTAIL, nodeIds.get(source)));
                    } else {
                        sourceNode.setName(nodeIds.get(source));
                    }
                    edgeStatement.setSourceNode(sourceNode);
                    // set target node or cluster
                    EdgeTarget edgeTarget = DotFactory.eINSTANCE.createEdgeTarget();
                    Node targetNode = DotFactory.eINSTANCE.createNode();
                    if (hierarchy && !target.getChildren().isEmpty()) {
                        targetNode.setName(target.getProperty(CLUSTER_DUMMY));
                        attributes.add(createAttribute(Attributes.LHEAD, nodeIds.get(target)));
                    } else {
                        targetNode.setName(nodeIds.get(target));
                    }
                    edgeTarget.setTargetnode(targetNode);
                    edgeStatement.getEdgeTargets().add(edgeTarget);

                    // add edge labels at head, tail, and middle position
                    setEdgeLabels(edge, attributes, vertical);
                    
                    if (transData.getProperty(USE_EDGE_IDS)) {
                        // add comment with edge identifier
                        String edgeID = getEdgeID(edge, transData);
                        attributes.add(createAttribute(Attributes.COMMENT, "\"" + edgeID + "\""));
                    }
                    
                    // include edge routing for full export, if there is one
                    if (!edge.getSections().isEmpty()) {
                        ElkEdgeSection edgeSection = edge.getSections().get(0);
                        
                        if (transformEdgeLayout && (edgeSection.getBendPoints().size() > 0
                                || edgeSection.getStartX() != 0 || edgeSection.getStartY() != 0
                                || edgeSection.getEndX() != 0 || edgeSection.getEndY() != 0)) {
                            
                            StringBuilder pos = new StringBuilder();
                            Iterator<KVector> pointIter = ElkUtil.createVectorChain(edgeSection).iterator();
                            while (pointIter.hasNext()) {
                                KVector point = pointIter.next();
                                ElkUtil.toAbsolute(point, edge.getContainingNode());
                                pos.append(point.x);
                                pos.append(",");
                                pos.append(point.y);
                                if (pointIter.hasNext()) {
                                    pos.append(" ");
                                }
                            }
                            attributes.add(createAttribute(Attributes.POS, "\"" + pos + "\""));
                        }
                    }
                    
                    statements.add(edgeStatement);
                }
            }
            if (hierarchy) {
                nodes.addAll(source.getChildren());
            }
        }
    }
    
    /**
     * Determines whether the given node is contained in the currently processed layout graph.
     * 
     * @param node a node
     * @param root the root node
     * @return true if the node is in the layout graph
     */
    private boolean isInsideGraph(final ElkNode node, final ElkNode root) {
        ElkNode n = node;
        do {
            if (n == root) {
                return true;
            }
            n = n.getParent();
        } while (n != null);
        return false;
    }

    /**
     * Sets attributes for the whole graph.
     * 
     * @param statements
     *            the statement list for adding attributes
     * @param parentNode
     *            the parent node
     * @param transData
     *            transformation data
     */
    protected void setGraphAttributes(final List<Statement> statements, final ElkNode parentNode,
            final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        // Implement in subclasses
    }

    /**
     * Sets general attributes for nodes.
     * 
     * @param statements
     *            the statement list for adding attributes
     * @return the list of node attributes
     */
    protected List<Attribute> setGeneralNodeAttributes(final List<Statement> statements) {
        AttributeStatement nodeAttrStatement = DotFactory.eINSTANCE.createAttributeStatement();
        nodeAttrStatement.setType(AttributeType.NODE);
        List<Attribute> nodeAttrs = nodeAttrStatement.getAttributes();
        statements.add(nodeAttrStatement);
        nodeAttrs.add(createAttribute(Attributes.SHAPE, "box"));
        nodeAttrs.add(createAttribute(Attributes.FIXEDSIZE, "true"));
        return nodeAttrs;
    }

    /**
     * Sets general attributes for edges.
     * 
     * @param statements
     *            the statement list for adding attributes
     * @return the list of edge attributes
     */
    protected List<Attribute> setGeneralEdgeAttributes(final List<Statement> statements) {
        AttributeStatement edgeAttrStatement = DotFactory.eINSTANCE.createAttributeStatement();
        edgeAttrStatement.setType(AttributeType.EDGE);
        List<Attribute> edgeAttrs = edgeAttrStatement.getAttributes();
        statements.add(edgeAttrStatement);
        edgeAttrs.add(createAttribute(Attributes.EDGEDIR, "none"));
        return edgeAttrs;
    }

    /**
     * Create an attribute with given name and value for the Dot graph.
     * 
     * @param name name of the attribute
     * @param value value of the attribute
     * @return instance of a Dot attribute
     */
    public static Attribute createAttribute(final String name, final String value) {
        Attribute attribute = DotFactory.eINSTANCE.createAttribute();
        attribute.setName(name);
        attribute.setValue(value);
        return attribute;
    }
    
    /**
     * Create an attribute with given name and integer value for the Dot graph.
     * 
     * @param name name of the attribute
     * @param value value of the attribute
     * @return instance of a Dot attribute
     */
    public static Attribute createAttribute(final String name, final int value) {
        Attribute attribute = DotFactory.eINSTANCE.createAttribute();
        attribute.setName(name);
        attribute.setValue("\"" + value + "\"");
        return attribute;
    }

    /**
     * Create an attribute with given name and float value for the Dot graph.
     * 
     * @param name name of the attribute
     * @param value value of the attribute
     * @return instance of a Dot attribute
     */
    public static Attribute createAttribute(final String name, final float value) {
        Attribute attribute = DotFactory.eINSTANCE.createAttribute();
        attribute.setName(name);
        attribute.setValue("\"" + value + "\"");
        return attribute;
    }

    /**
     * Create an attribute with given name and float value for the Dot graph.
     * 
     * @param name name of the attribute
     * @param value value of the attribute
     * @return instance of a Dot attribute
     */
    public static Attribute createAttribute(final String name, final double value) {
        Attribute attribute = DotFactory.eINSTANCE.createAttribute();
        attribute.setName(name);
        attribute.setValue("\"" + value + "\"");
        return attribute;
    }

    /**
     * Set edge labels for the given edge.
     * 
     * @param elkedge edge whose labels shall be set
     * @param attributes edge attribute list to which the labels are added
     * @param isVertical indicates whether vertical layout direction is active
     */
    protected void setEdgeLabels(final ElkEdge elkedge, final List<Attribute> attributes,
            final boolean isVertical) {
        
        if (elkedge.getLabels().isEmpty()) {
            return;
        }
        // as Graphviz only supports positioning of one label per label placement, all labels
        // are stacked to one big label as workaround
        StringBuilder midLabel = new StringBuilder(), headLabel = new StringBuilder(),
                tailLabel = new StringBuilder();
        String fontName = null;
        int fontSize = 0;
        boolean isCenterFontName = false, isCenterFontSize = false;
        for (ElkLabel label : elkedge.getLabels()) {
            StringBuilder buffer = midLabel;
            EdgeLabelPlacement placement = label.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT);
            boolean takeFontName = false, takeFontSize = false;
            switch (placement) {
            case HEAD:
                takeFontName = fontName == null;
                takeFontSize = fontSize <= 0;
                buffer = headLabel;
                break;
            case TAIL:
                takeFontName = fontName == null;
                takeFontSize = fontSize <= 0;
                buffer = tailLabel;
                break;
            default: // CENTER
                takeFontName = fontName == null || !isCenterFontName;
                isCenterFontName = true;
                takeFontSize = fontSize <= 0 || !isCenterFontSize;
                isCenterFontSize = true;
                break;
            }
            if (buffer.length() > 0) {
                buffer.append("\n");
            }
            buffer.append(label.getText());
            if (takeFontName) {
                fontName = label.getProperty(CoreOptions.FONT_NAME);
            }
            if (takeFontSize) {
                fontSize = label.getProperty(CoreOptions.FONT_SIZE);
                // increase the font size to let Graphviz prepare more room for the label
                if (fontSize > 0) {
                    fontSize *= FONT_SIZE_MULT;
                }
            }
        }

        // set mid label: if empty, it is filled with a dummy string to avoid
        // edge overlapping
        if (midLabel.length() > 0) {
            double labelSpacing = elkedge.getProperty(CoreOptions.SPACING_EDGE_LABEL);
            if (labelSpacing < 1) {
                labelSpacing = 0;
            }
            int charsToAdd = (int) labelSpacing - 1;
            for (int i = 0; i < charsToAdd; i++) {
                midLabel.append(isVertical ? "O" : "\nO");
            }
            attributes.add(createAttribute(Attributes.LABEL, createString(midLabel.toString())));
        }
        // set head label
        if (headLabel.length() > 0) {
            attributes.add(createAttribute(Attributes.HEADLABEL, createString(headLabel.toString())));
        }
        // set tail label
        if (tailLabel.length() > 0) {
            attributes.add(createAttribute(Attributes.TAILLABEL, createString(tailLabel.toString())));
        }
        // set font name
        if (fontName != null && fontName.length() > 0) {
            attributes.add(createAttribute(Attributes.FONTNAME, "\"" + fontName + "\""));
        }
        // set font size
        if (fontSize > 0) {
            attributes.add(createAttribute(Attributes.FONTSIZE, fontSize));
        }
    }

    /** first character that is not replaced by underscore. */
    private static final char MIN_OUT_CHAR = 32;
    /** last character that is not replaced by underscore. */
    private static final char MAX_OUT_CHAR = 126;

    /**
     * Creates a properly parseable string by adding the escape character '\\' wherever needed and
     * replacing illegal characters.
     * 
     * @param label a label string from the KGraph structure
     * @return string to be used in the Graphviz model
     */
    private static String createString(final String label) {
        StringBuilder escapeBuffer = new StringBuilder(label.length() + 2);
        // prefix the label with an underscore to prevent it from being equal to
        // a keyword
        escapeBuffer.append("\"");
        for (int i = 0; i < label.length(); i++) {
            char c = label.charAt(i);
            if (c == '\"' || c == '\\' || c > MAX_OUT_CHAR) {
                escapeBuffer.append('_');
            } else if (c == '\n') {
                escapeBuffer.append("\\n");
            } else if (c >= MIN_OUT_CHAR) {
                escapeBuffer.append(c);
            }
        }
        escapeBuffer.append('\"');
        return escapeBuffer.toString();
    }

    /** node types used for identifier generation. */
    private enum NodeType {
        NODE, CLUSTER, DUMMY;
    }
    
    /**
     * Creates a unique identifier for the given node.
     * 
     * @param node node for which an identifier shall be created
     * @param cluster whether a cluster id should be created
     * @param transData transformation data
     * @return a unique string used to identify the node
     */
    private String getNodeID(final ElkNode node, final NodeType type,
            final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        int id = transData.getProperty(NEXT_NODE_ID);
        transData.setProperty(NEXT_NODE_ID, id + 1);
        String idstring = null;
        switch (type) {
        case NODE:
            idstring = "node" + id;
            transData.getProperty(GRAPH_ELEMS).put(idstring, node);
            break;
        case CLUSTER:
            idstring = "cluster" + id;
            transData.getProperty(GRAPH_ELEMS).put(idstring, node);
            break;
        case DUMMY:
            idstring = "dummy" + id;
            node.setProperty(CLUSTER_DUMMY, idstring);
            break;
        }
        return idstring;
    }

    /**
     * Creates a unique identifier for the given edge.
     * 
     * @param edge edge for which an identifier shall be created
     * @param transData transformation data
     * @return a unique string used to identify the edge
     */
    private String getEdgeID(final ElkEdge edge, final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        int id = transData.getProperty(NEXT_EDGE_ID);
        transData.setProperty(NEXT_EDGE_ID, id + 1);
        String idstring = "edge" + id;
        transData.getProperty(GRAPH_ELEMS).put(idstring, edge);
        return idstring;
    }

    /**
     * Applies layout information from a Graphviz model to the original graph.
     * Note that GraphViz uses cubic B-Splines for edge routing, some generalization of Bezier curves.
     * Edge offsets are given separately, since on some platforms edges seem to have a different
     * reference point than nodes.
     * 
     * @param parentNode parent node of the original graph
     * @param statements list of statements to process
     * @param baseOffset offset to be added to edge and subgraph coordinates
     * @param transData transformation data
     */
    private void applyLayout(final ElkNode parentNode, final List<Statement> statements, final KVector baseOffset,
            final ElkPadding outerPadding, final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        // process attributes: determine bounding box of the parent node
        ElkPadding padding = outerPadding;
        KVector nodeOffset = new KVector();
        attr_loop: for (Statement statement : statements) {
            if (statement instanceof AttributeStatement) {
                AttributeStatement attributeStatement = (AttributeStatement) statement;
                if (attributeStatement.getType() == AttributeType.GRAPH) {
                    for (Attribute attribute : attributeStatement.getAttributes()) {
                        if (attribute.getName().equals(Attributes.BOUNDINGBOX)) {
                            try {
                                StringTokenizer tokenizer = new StringTokenizer(attribute.getValue(),
                                        ATTRIBUTE_DELIM);
                                double leftx = Double.parseDouble(tokenizer.nextToken());
                                double bottomy = Double.parseDouble(tokenizer.nextToken());
                                double rightx = Double.parseDouble(tokenizer.nextToken());
                                double topy = Double.parseDouble(tokenizer.nextToken());
                                // set parent node attributes
                                double width = rightx - leftx;
                                double height = bottomy - topy;
                                
                                if (parentNode == transData.getSourceGraph()) {
                                    // the root node, representing the drawing frame
                                    width += padding.getHorizontal();
                                    height += padding.getVertical();
                                    baseOffset.add(-leftx, -topy);
                                    nodeOffset.add(-leftx, -topy);
                                } else {
                                    // a child or descendant of the root node 
                                    parentNode.setX(baseOffset.x + leftx + padding.getLeft());
                                    parentNode.setY(baseOffset.y + topy + padding.getTop());
                                    
                                    // since dot uses a 'compound' layout instead of laying out a hierarchical graph, 
                                    // no padding is supported for children
                                    padding = new ElkPadding();
                                    // ... and the children's offset must be somewhat adjusted
                                    nodeOffset.x = -(baseOffset.x + leftx);
                                    nodeOffset.y = -(baseOffset.y + topy);
                                }
                                ElkUtil.resizeNode(parentNode, width, height, false, true);
                                parentNode.setProperty(CoreOptions.NODE_SIZE_CONSTRAINTS, SizeConstraint.fixed());
                                break attr_loop;
                            } catch (NumberFormatException exception) {
                                // ignore exception
                            } catch (NoSuchElementException exception) {
                                // ignore exception
                            }
                        }
                    }
                }
            }
        }
        
        // process nodes and subgraphs to collect all offset data
        for (Statement statement : statements) {
            if (statement instanceof NodeStatement) {
                applyNodeLayout((NodeStatement) statement, nodeOffset, padding, transData);
            } else if (statement instanceof Subgraph) {
                Subgraph subgraph = (Subgraph) statement;
                ElkNode elknode = (ElkNode) transData.getProperty(GRAPH_ELEMS).get(subgraph.getName());
                applyLayout(elknode, subgraph.getStatements(), baseOffset, padding, transData);
                elknode.setX(elknode.getX() + nodeOffset.x);
                elknode.setY(elknode.getY() + nodeOffset.y);
            }
        }
    }
    
    /**
     * Set the position of a node.
     * 
     * @param nodeStatement a node statement
     * @param offset the offset to be added to node coordinates
     * @param spacing additional border spacing
     * @param transData transformation data
     */
    private void applyNodeLayout(final NodeStatement nodeStatement, final KVector offset,
            final ElkPadding padding, final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        ElkNode elknode = (ElkNode) transData.getProperty(GRAPH_ELEMS).get(nodeStatement.getNode().getName());
        if (elknode == null) {
            return;
        }
        
        double xpos = 0.0;
        double ypos = 0.0;
        double width = 0.0;
        double height = 0.0f;
        for (Attribute attribute : nodeStatement.getAttributes()) {
            try {
                if (attribute.getName().equals(Attributes.POS)) {
                    KVector pos = new KVector();
                    pos.parse((attribute.getValue()));
                    xpos = pos.x + offset.x + padding.getLeft();
                    ypos = pos.y + offset.y + padding.getTop();
                } else if (attribute.getName().equals(Attributes.WIDTH)) {
                    StringTokenizer tokenizer = new StringTokenizer(attribute.getValue(), ATTRIBUTE_DELIM);
                    width = Double.parseDouble(tokenizer.nextToken()) * DPI;
                } else if (attribute.getName().equals(Attributes.HEIGHT)) {
                    StringTokenizer tokenizer = new StringTokenizer(attribute.getValue(), ATTRIBUTE_DELIM);
                    height = Double.parseDouble(tokenizer.nextToken()) * DPI;
                }
            } catch (NumberFormatException exception) {
                // ignore exception
            } catch (IllegalArgumentException exception) {
                // ignore exception
            } catch (NoSuchElementException exception) {
                // ignore exception
            }
        }
        elknode.setX(xpos - width / 2);
        elknode.setY(ypos - height / 2);
    }
    
    /**
     * Sets the bend points of an edge.
     * 
     * @param edgeStatement an edge statement
     * @param edgeOffset offset to be added to edge coordinates
     * @param transData transformation data
     */
    private void applyEdgeLayout(final EdgeStatement edgeStatement, final KVector edgeOffset,
            final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        Map<String, String> attributeMap = createAttributeMap(edgeStatement.getAttributes());
        ElkEdge elkedge = (ElkEdge) transData.getProperty(GRAPH_ELEMS).get(attributeMap.get(Attributes.COMMENT));
        if (elkedge == null) {
            return;
        }
        
        ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(elkedge, true, true);
        String posString = attributeMap.get(Attributes.POS);
        if (posString == null) {
            posString = "";
        }
        
        ElkNode referenceNode = elkedge.getContainingNode();
        KVector reference = new KVector();
        while (referenceNode != null && referenceNode != transData.getSourceGraph()) {
            // Note: do _not_ apply padding here,
            // it is already included in the the 'edgeOffset'
            reference.x += referenceNode.getX();
            reference.y += referenceNode.getY();
            referenceNode = referenceNode.getParent();
        }
        KVector offset = edgeOffset.clone().sub(reference);

        // parse the list of spline control coordinates
        List<KVectorChain> splines = new LinkedList<KVectorChain>();
        Pair<KVector, KVector> endpoints = parseSplinePoints(posString, splines, offset);

        KVector sourcePoint = endpoints.getFirst();
        KVector targetPoint = endpoints.getSecond();
        if (!splines.isEmpty()) {
            // the first point in the list is the start point, if no arrowhead is given
            if (sourcePoint == null) {
                List<KVector> points = splines.get(0);
                if (!points.isEmpty()) {
                    sourcePoint = points.remove(0);
                } else {
                    ElkNode sourceNode = ElkGraphUtil.connectableShapeToNode(elkedge.getSources().get(0));
                    sourcePoint = new KVector();
                    sourcePoint.x = sourceNode.getX() + sourceNode.getWidth() / 2;
                    sourcePoint.y = sourceNode.getY() + sourceNode.getHeight() / 2;
                }
            }
            
            // the last point in the list is the end point, if no arrowhead is given
            if (targetPoint == null) {
                List<KVector> points = splines.get(splines.size() - 1);
                if (!points.isEmpty()) {
                    targetPoint = points.remove(points.size() - 1);
                } else {
                    ElkNode targetNode = ElkGraphUtil.connectableShapeToNode(elkedge.getTargets().get(0));
                    targetPoint = new KVector();
                    targetPoint.x = targetNode.getX() + targetNode.getWidth() / 2;
                    targetPoint.y = targetNode.getY() + targetNode.getHeight() / 2;
                }
            }
            
            // add all other control points to the edge
            for (KVectorChain points : splines) {
                for (KVector point : points) {
                    ElkGraphUtil.createBendPoint(edgeSection, point.x, point.y);
                }
            }
            edgeSection.setStartLocation(sourcePoint.x, sourcePoint.y);
            edgeSection.setEndLocation(targetPoint.x, targetPoint.y);
            
            // set source and target port positions accordingly
            boolean adaptPortPositions = transData.getProperty(ADAPT_PORT_POSITIONS);
            if (adaptPortPositions
                    && (elkedge.getSources().get(0) instanceof ElkPort
                            || elkedge.getTargets().get(0) instanceof ElkPort)) {
                
                referenceNode = elkedge.getContainingNode();
                
                if (elkedge.getSources().get(0) instanceof ElkPort) {
                    ElkPort sourcePort = (ElkPort) elkedge.getSources().get(0);
                    ElkNode sourceNode = sourcePort.getParent();
                    
                    ElkUtil.toAbsolute(sourcePoint, referenceNode);
                    ElkUtil.toRelative(sourcePoint, sourceNode);
                    
                    sourcePort.setX((float) sourcePoint.x - sourcePort.getWidth() / 2);
                    sourcePort.setY((float) sourcePoint.y - sourcePort.getHeight() / 2);
                }
                
                if (elkedge.getTargets().get(0) instanceof ElkPort) {
                    ElkPort targetPort = (ElkPort) elkedge.getTargets().get(0);
                    ElkNode targetNode = targetPort.getParent();
                    
                    ElkUtil.toAbsolute(targetPoint, referenceNode);
                    ElkUtil.toRelative(targetPoint, targetNode);
                    
                    targetPort.setX(targetPoint.x - targetPort.getWidth() / 2);
                    targetPort.setY(targetPoint.y - targetPort.getHeight() / 2);
                }
            }
        }
        
        if (transData.getProperty(USE_SPLINES)) {
            elkedge.setProperty(CoreOptions.EDGE_ROUTING, EdgeRouting.SPLINES);
        }

        // process the edge labels
        String labelPos = attributeMap.get(Attributes.LABELPOS);
        if (labelPos != null) {
            applyEdgeLabelPos(elkedge, labelPos, EdgeLabelPlacement.CENTER, offset);
        }
        labelPos = attributeMap.get(Attributes.HEADLP);
        if (labelPos != null) {
            applyEdgeLabelPos(elkedge, labelPos, EdgeLabelPlacement.HEAD, offset);
        }
        labelPos = attributeMap.get(Attributes.TAILLP);
        if (labelPos != null) {
            applyEdgeLabelPos(elkedge, labelPos, EdgeLabelPlacement.TAIL, offset);
        }
    }

    /**
     * Applies the edge label positions for the given edge.
     * 
     * @param elkedge edge for which labels are processed
     * @param posString string with label position
     * @param placement label placement to choose
     * @param offsetx x offset added to positions
     * @param offsety y offset added to positions
     */
    private void applyEdgeLabelPos(final ElkEdge elkedge, final String posString,
            final EdgeLabelPlacement placement, final KVector offset) {
        
        double combinedWidth = 0.0;
        double combinedHeight = 0.0;
        for (ElkLabel label : elkedge.getLabels()) {
            EdgeLabelPlacement elp = label.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT);
            if (elp == placement) {
                combinedWidth = Math.max(combinedWidth, label.getWidth());
                combinedHeight += label.getHeight();
            }
        }
        try {
            KVector pos = new KVector();
            pos.parse(posString);
            
            double xpos = pos.x - combinedWidth / 2 + offset.x;
            double ypos = pos.y - combinedHeight / 2 + offset.y;
            
            for (ElkLabel label : elkedge.getLabels()) {
                EdgeLabelPlacement elp = label.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT);
                if (elp == placement) {
                    double xoffset = (combinedWidth - label.getWidth()) / 2;
                    label.setX(xpos + xoffset);
                    label.setY(ypos);
                    ypos += label.getHeight();
                }
            }
        } catch (IllegalArgumentException exception) {
            // ignore exception
        }
    }

    /**
     * Converts a list of attributes to a map of attribute names to their values.
     * 
     * @param attributes attribute list
     * @return a hash map that contains all given attributes
     */
    private static Map<String, String> createAttributeMap(final List<Attribute> attributes) {
        Map<String, String> attributeMap = new HashMap<String, String>(attributes.size());
        for (Attribute attribute : attributes) {
            attributeMap.put(attribute.getName(), attribute.getValue());
        }
        return attributeMap;
    }

    /**
     * Puts the points of a position string into a list of splines.
     * 
     * @param posString string with spline points
     * @param splines list of splines
     * @param offset offset to add to coordinates
     * @return the source and the target point, if specified by the position string
     */
    private static Pair<KVector, KVector> parseSplinePoints(final String posString,
            final List<KVectorChain> splines, final KVector offset) {
        KVector sourcePoint = null, targetPoint = null;
        StringTokenizer splinesTokenizer = new StringTokenizer(posString, "\";");
        while (splinesTokenizer.hasMoreTokens()) {
            KVectorChain pointList = new KVectorChain();
            StringTokenizer posTokenizer = new StringTokenizer(splinesTokenizer.nextToken(), " \t");
            while (posTokenizer.hasMoreTokens()) {
                String token = posTokenizer.nextToken();
                try {
                    if (token.startsWith("s")) {
                        if (sourcePoint == null) {
                            sourcePoint = new KVector();
                            int commaIndex = token.indexOf(',');
                            sourcePoint.parse(token.substring(commaIndex + 1));
                            sourcePoint.add(offset);
                        }
                    } else if (token.startsWith("e")) {
                        if (targetPoint == null) {
                            targetPoint = new KVector();
                            int commaIndex = token.indexOf(',');
                            targetPoint.parse(token.substring(commaIndex + 1));
                            targetPoint.add(offset);
                        }
                    } else {
                        KVector point = new KVector();
                        point.parse(token);
                        pointList.add(point.add(offset));
                    }
                } catch (IllegalArgumentException exception) {
                    // ignore exception
                }
            }
            splines.add(pointList);
        }
        return new Pair<KVector, KVector>(sourcePoint, targetPoint);
    }
    
}
