/**
 * Copyright (c) 2016 Kiel University and others.
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
 */
package org.eclipse.elk.graph.util;

import java.util.Map;

import org.eclipse.elk.graph.*;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

import org.eclipse.emf.common.notify.Adapter;
import org.eclipse.emf.common.notify.Notifier;

import org.eclipse.emf.common.notify.impl.AdapterFactoryImpl;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * The <b>Adapter Factory</b> for the model.
 * It provides an adapter <code>createXXX</code> method for each class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.graph.ElkGraphPackage
 * @generated
 */
public class ElkGraphAdapterFactory extends AdapterFactoryImpl {
    /**
     * The cached model package.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected static ElkGraphPackage modelPackage;

    /**
     * Creates an instance of the adapter factory.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkGraphAdapterFactory() {
        if (modelPackage == null) {
            modelPackage = ElkGraphPackage.eINSTANCE;
        }
    }

    /**
     * Returns whether this factory is applicable for the type of the object.
     * <!-- begin-user-doc -->
     * This implementation returns <code>true</code> if the object is either the model's package or is an instance object of the model.
     * <!-- end-user-doc -->
     * @return whether this factory is applicable for the type of the object.
     * @generated
     */
    @Override
    public boolean isFactoryForType(Object object) {
        if (object == modelPackage) {
            return true;
        }
        if (object instanceof EObject) {
            return ((EObject)object).eClass().getEPackage() == modelPackage;
        }
        return false;
    }

    /**
     * The switch that delegates to the <code>createXXX</code> methods.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected ElkGraphSwitch<Adapter> modelSwitch =
        new ElkGraphSwitch<Adapter>() {
            @Override
            public Adapter caseIPropertyHolder(IPropertyHolder object) {
                return createIPropertyHolderAdapter();
            }
            @Override
            public Adapter caseEMapPropertyHolder(EMapPropertyHolder object) {
                return createEMapPropertyHolderAdapter();
            }
            @Override
            public Adapter caseElkGraphElement(ElkGraphElement object) {
                return createElkGraphElementAdapter();
            }
            @Override
            public Adapter caseElkShape(ElkShape object) {
                return createElkShapeAdapter();
            }
            @Override
            public Adapter caseElkLabel(ElkLabel object) {
                return createElkLabelAdapter();
            }
            @Override
            public Adapter caseElkConnectableShape(ElkConnectableShape object) {
                return createElkConnectableShapeAdapter();
            }
            @Override
            public Adapter caseElkNode(ElkNode object) {
                return createElkNodeAdapter();
            }
            @Override
            public Adapter caseElkPort(ElkPort object) {
                return createElkPortAdapter();
            }
            @Override
            public Adapter caseElkEdge(ElkEdge object) {
                return createElkEdgeAdapter();
            }
            @Override
            public Adapter caseElkBendPoint(ElkBendPoint object) {
                return createElkBendPointAdapter();
            }
            @Override
            public Adapter caseElkEdgeSection(ElkEdgeSection object) {
                return createElkEdgeSectionAdapter();
            }
            @Override
            public Adapter caseElkPropertyToValueMapEntry(Map.Entry<IProperty<?>, Object> object) {
                return createElkPropertyToValueMapEntryAdapter();
            }
            @Override
            public Adapter defaultCase(EObject object) {
                return createEObjectAdapter();
            }
        };

    /**
     * Creates an adapter for the <code>target</code>.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param target the object to adapt.
     * @return the adapter for the <code>target</code>.
     * @generated
     */
    @Override
    public Adapter createAdapter(Notifier target) {
        return modelSwitch.doSwitch((EObject)target);
    }


    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.properties.IPropertyHolder <em>IProperty Holder</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.properties.IPropertyHolder
     * @generated
     */
    public Adapter createIPropertyHolderAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.EMapPropertyHolder <em>EMap Property Holder</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.EMapPropertyHolder
     * @generated
     */
    public Adapter createEMapPropertyHolderAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkGraphElement <em>Element</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkGraphElement
     * @generated
     */
    public Adapter createElkGraphElementAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkShape <em>Elk Shape</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkShape
     * @generated
     */
    public Adapter createElkShapeAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkLabel <em>Elk Label</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkLabel
     * @generated
     */
    public Adapter createElkLabelAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkConnectableShape <em>Elk Connectable Shape</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkConnectableShape
     * @generated
     */
    public Adapter createElkConnectableShapeAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkNode <em>Elk Node</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkNode
     * @generated
     */
    public Adapter createElkNodeAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkPort <em>Elk Port</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkPort
     * @generated
     */
    public Adapter createElkPortAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkEdge <em>Elk Edge</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkEdge
     * @generated
     */
    public Adapter createElkEdgeAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkBendPoint <em>Elk Bend Point</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkBendPoint
     * @generated
     */
    public Adapter createElkBendPointAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link org.eclipse.elk.graph.ElkEdgeSection <em>Elk Edge Section</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see org.eclipse.elk.graph.ElkEdgeSection
     * @generated
     */
    public Adapter createElkEdgeSectionAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for an object of class '{@link java.util.Map.Entry <em>Elk Property To Value Map Entry</em>}'.
     * <!-- begin-user-doc -->
     * This default implementation returns null so that we can easily ignore cases;
     * it's useful to ignore a case when inheritance will catch all the cases anyway.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @see java.util.Map.Entry
     * @generated
     */
    public Adapter createElkPropertyToValueMapEntryAdapter() {
        return null;
    }

    /**
     * Creates a new adapter for the default case.
     * <!-- begin-user-doc -->
     * This default implementation returns null.
     * <!-- end-user-doc -->
     * @return the new adapter.
     * @generated
     */
    public Adapter createEObjectAdapter() {
        return null;
    }

} //ElkGraphAdapterFactory
